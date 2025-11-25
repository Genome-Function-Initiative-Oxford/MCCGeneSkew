import argparse
from glob import glob
import os
import pandas as pd
import pysam
from collections import namedtuple
from scipy.stats import binomtest

ReadEvidence = namedtuple("ReadEvidence", "ref alt")


class CombinedReadEvidence(object):
    def __init__(
        self, up, down, both, up_no_indel, down_no_indel
    ):
        self.up = up
        self.down = down
        self.both = both
        self.up_no_indel = up_no_indel
        self.down_no_indel = down_no_indel

    @property
    def ref(self):
        return self.both.ref

    @property
    def alt(self):
        return self.both.alt

    def binomtest(self):
        if self.ref + self.alt >= 2:
            return binomtest(self.ref // 2, (self.ref + self.alt) // 2)
        else:
            return binomtest(1, 2)


def get_hets(peak, variant_file_lookup):
    hets = set()
    for donor, vf in variant_file_lookup.items():
        hets.update({
            ":".join((v.chrom, str(v.pos), v.ref, v.alts[0])) for v in vf.fetch(
                contig=peak.Chromosome,
                start=peak.Start,
                end=peak.End
            ) if set(v.samples[donor]["GT"]) == {0, 1} and "PASS" in v.filter
        })
    return sorted(hets, key=lambda x: int(x.split(":")[1]))


def get_ref_allele_count(variant, donor, phase_index, key="A"):
    if variant.samples[donor]["GT"][phase_index] == 0:
        return variant.info.get(key + "R", 0)
    else:
        return variant.info.get(key + "A", 0)


def get_alt_allele_count(variant, donor, phase_index, key="A"):
    if variant.samples[donor]["GT"][phase_index] == 1:
        return variant.info.get(key + "R", 0)
    else:
        return variant.info.get(key + "A", 0)
    

def filter_variants(variants):
    current_pos = 0
    filtered_variants = []
    for v in variants:
        if v.start >= current_pos + 40:
            filtered_variants.append(v)
            current_pos = v.stop
    return filtered_variants

    
def get_annotations(het, donor, variant_file, start, end, key):
    chrom = het.split(":")[0]
    pos = int(het.split(":")[1])
    ref = het.split(":")[2]
    alt = het.split(":")[3]

    variants = list(variant_file.fetch(chrom, start, end))

    has_variant = len([v for v in variants if v.pos == pos and v.ref == ref and alt in v.alts]) > 0
    if not has_variant:
        return None

    main_variant = [v for v in variants if v.pos == pos and v.ref == ref and alt in v.alts][0]
    if not set(main_variant.samples[donor]["GT"]) == {0, 1}:
        return None

    # Compute phase of alternative
    phase_index = 0 if main_variant.samples[donor]["GT"][0] == 0 else 1

    phase_set = main_variant.samples[donor]["PS"]
    phase_quality = main_variant.samples[donor].get("PQ", 0)

    min_phase_qual = 10
    if phase_quality >= min_phase_qual:
        phased_variants = [
            v for v in variants if v.samples[donor]["PS"] == phase_set 
            and set(v.samples[donor]["GT"]) == {0, 1} and
            v.samples[donor].get("PQ", 0) >= min_phase_qual
        ]
    else:
        phased_variants = [main_variant]

    aligned_up = [v for v in phased_variants if v.samples[donor]["GT"] == (1, 0)]
    aligned_down = [v for v in phased_variants if v.samples[donor]["GT"] == (0, 1)]

    aligned_up_no_indel = [v for v in aligned_up if len(v.ref) == 1 and all((len(alt[0]) == 1 for alt in v.alts))]
    aligned_down_no_indel = [v for v in aligned_down if len(v.ref) == 1 and all((len(alt[0]) == 1 for alt in v.alts))]

    return CombinedReadEvidence(
        up=ReadEvidence(
            sum([get_ref_allele_count(v, donor, phase_index, key=key) for v in filter_variants(aligned_up)]),
            sum([get_alt_allele_count(v, donor, phase_index, key=key) for v in filter_variants(aligned_up)])
        ),
        down=ReadEvidence(
            sum([get_ref_allele_count(v, donor, phase_index, key=key) for v in filter_variants(aligned_down)]),
            sum([get_alt_allele_count(v, donor, phase_index, key=key) for v in filter_variants(aligned_down)])
        ),
        up_no_indel=ReadEvidence(
            sum([get_ref_allele_count(v, donor, phase_index, key=key) for v in filter_variants(aligned_up_no_indel)]),
            sum([get_alt_allele_count(v, donor, phase_index, key=key) for v in filter_variants(aligned_up_no_indel)])
        ),
        down_no_indel=ReadEvidence(
            sum([get_ref_allele_count(v, donor, phase_index, key=key) for v in filter_variants(aligned_down_no_indel)]),
            sum([get_alt_allele_count(v, donor, phase_index, key=key) for v in filter_variants(aligned_down_no_indel)])
        ),
        both=ReadEvidence(
            sum([get_ref_allele_count(v, donor, phase_index, key=key) for v in filter_variants(phased_variants)]),
            sum([get_alt_allele_count(v, donor, phase_index, key=key) for v in filter_variants(phased_variants)])
        )
    )


def get_bam_data(het, variant_file_lookup, peak_start, peak_end):  
    atac = {}
    ctcf = {}
    k27 = {}
    k41 = {}
    k43 = {}
    
    for donor, vcf_file in variant_file_lookup.items():
        annotations_atac = get_annotations(
            het, donor, vcf_file, peak_start, peak_end, "A"
        )
        
        if annotations_atac is None:
            continue
        
        atac[donor] = annotations_atac
        ctcf[donor] = get_annotations(het, donor, vcf_file, peak_start, peak_end, "C")

        no_evidence = CombinedReadEvidence(
            ReadEvidence(0, 0), ReadEvidence(0, 0), ReadEvidence(0, 0), ReadEvidence(0, 0), ReadEvidence(0, 0)
        )
        
        mid_peak = (peak_start + peak_end) // 2
        broad_start = min(peak_start, mid_peak - 1500)
        broad_end = max(peak_end, mid_peak + 1500)

        k27[donor] = get_annotations(
            het, donor, vcf_file, broad_start, broad_end, "K27"
        ) if "K27A" in vcf_file.header.info.keys() else no_evidence
        
        k41[donor] = get_annotations(
            het, donor, vcf_file, broad_start, broad_end, "K41"
        ) if "K41A" in vcf_file.header.info.keys() else no_evidence
         
        k43[donor] = get_annotations(
            het, donor, vcf_file, broad_start, broad_end, "K43"
        ) if "K43A" in vcf_file.header.info.keys() else no_evidence
    
    return atac, ctcf, k27, k41, k43


def sum_ref(count_lookup):
    return sum([x.ref for x in count_lookup.values()])


def sum_alt(count_lookup):
    return sum([x.alt for x in count_lookup.values()])


if __name__=="__main__":
    parser = argparse.ArgumentParser(
        prog='CollateSkew',
        description='Takes Donor + bam file patterns. Builds config for analysis',
    )
    parser.add_argument(
        "--regions-file",
        default="/project/Wellcome_Discovery/lhentges/atac_LoTronOpV/union_with_peak_score_plus_pvals.bed"
    )
    parser.add_argument(
        "--regions-black-list",
        default="/project/Wellcome_Discovery/esanders/exclusion_regions/ENCFF356LFX.bed.gz"
    )
    parser.add_argument(
        "--skew-file-pattern",
        default="/project/Wellcome_Discovery/esanders/00_atac_skew/00_ATAC_counts/*.BAM_counts.vcf.gz"
    )
    parser.add_argument("--chrom")
    parser.add_argument("--output")

    args=parser.parse_args()

    regions = pd.read_csv(args.regions_file, sep="\t")
    regions = regions[regions.Chromosome == args.chrom]

    skew_files = glob(args.skew_file_pattern)
    print("Found {} files".format(len(skew_files)))
    
    variant_file_lookup = {}
    for f in skew_files:
        vf = pysam.VariantFile(f)
        sample_name = list(vf.header.samples)[0]
        variant_file_lookup[sample_name] = vf


    output = []
    for _, peak in regions.iterrows():
        hets = get_hets(peak, variant_file_lookup)

        for het in hets:
            atac, ctcf, k27, k41, k43 = get_bam_data(het, variant_file_lookup, peak.Start, peak.End)

            atac_ref = sum_ref(atac)
            atac_alt = sum_alt(atac)
            
            if atac_ref + atac_alt >= 2:
                test_stastic = binomtest(atac_ref // 2, (atac_ref + atac_alt) // 2)

                ref_high = atac_ref > atac_alt

                donor_confidence_intervals = {
                    donor: value.binomtest().proportion_ci(0.7) for donor, value in atac.items()
                }
                low_donors = [donor for donor, value in donor_confidence_intervals.items() if value.high < 0.5]
                unclear_donors = [
                    donor for donor, value in donor_confidence_intervals.items() if value.low <= 0.5 <= value.high]
                high_donors = [donor for donor, value in donor_confidence_intervals.items() if value.low > 0.5]

                output_item = {
                    "Chromsome": peak.Chromosome,
                    "Start": peak.Start,
                    "End": peak.End,
                    "het": het
                }
                
                output_item["ATAC_ref"] = atac_ref
                output_item["ATAC_alt"] = atac_alt
                
                output_item["CTCF_ref"] = sum_ref(ctcf)
                output_item["CTCF_alt"] = sum_alt(ctcf)

                output_item["H3K27ac_ref"] = sum_ref(k27)
                output_item["H3K27ac_alt"] = sum_alt(k27)

                output_item["H3K4me1_ref"] = sum_ref(k41)
                output_item["H3K4me1_alt"] = sum_alt(k41)

                output_item["H3K4me3_ref"] = sum_ref(k43)
                output_item["H3K4me3_alt"] = sum_alt(k43)

                output_item["ATAC_ref_up"] = sum([x.up.ref for x in atac.values()])
                output_item["ATAC_alt_up"] = sum([x.up.alt for x in atac.values()])
                
                output_item["ATAC_ref_down"] = sum([x.down.ref for x in atac.values()])
                output_item["ATAC_alt_down"] = sum([x.down.alt for x in atac.values()])

                output_item["CTCF_ref_up"] = sum([x.up.ref for x in ctcf.values()])
                output_item["CTCF_alt_up"] = sum([x.up.alt for x in ctcf.values()])
                
                output_item["CTCF_ref_down"] = sum([x.down.ref for x in ctcf.values()])
                output_item["CTCF_alt_down"] = sum([x.down.alt for x in ctcf.values()])

                output_item["H3K27ac_ref_down"] = sum([x.down.ref for x in k27.values()])
                output_item["H3K27ac_alt_down"] = sum([x.down.alt for x in k27.values()])

                output_item["H3K27ac_ref_up"] = sum([x.up.ref for x in k27.values()])
                output_item["H3K27ac_alt_up"] = sum([x.up.alt for x in k27.values()])

                output_item["ATAC_ref_up_no_indel"] = sum([x.up_no_indel.ref for x in atac.values()])
                output_item["ATAC_alt_up_no_indel"] = sum([x.up_no_indel.alt for x in atac.values()])

                output_item["ATAC_ref_down_no_indel"] = sum([x.down_no_indel.ref for x in atac.values()])
                output_item["ATAC_alt_down_no_indel"] = sum([x.down_no_indel.alt for x in atac.values()])

                output_item["For"] = ",".join(sorted(high_donors)) if ref_high else  ",".join(sorted(low_donors)) 
                output_item["Against"] = ",".join(sorted(low_donors)) if  ref_high else  ",".join(sorted(high_donors))
                output_item["Unclear"] = ",".join(sorted(unclear_donors))
                
                output.append(output_item)

    pd.DataFrame(output).to_csv(args.output, sep="\t", index=False)
