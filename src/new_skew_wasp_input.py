import argparse
import pysam
import pandas as pd
import json
import sys
from io import StringIO
import numpy as np


def get_variants(chrom, start, end, vf):
    return [
        v for v in vf.fetch(
            contig=chrom,
            start=start,
            end=end
        )
    ]


def get_haplotype_counts(allele_index, het_variants, donor, experiment_key):    
    return [
        v.info.get(experiment_key + "R", 0)
        if v.samples[donor]["GT"][0] == allele_index
        else v.info.get(experiment_key + "A", 0)
        for v in het_variants
    ]

def get_haplotypes_counts(het_variants, donor, experiment_key, has_allele_specific_alignment):
    
    haplotype_counts_0 = get_haplotype_counts(0, het_variants, donor, experiment_key)
    haplotype_counts_1 = get_haplotype_counts(1, het_variants, donor, experiment_key)
    if not has_allele_specific_alignment:
        # No personalised genomes for this data-type
        return haplotype_counts_0, haplotype_counts_1
    
    haplotype_counts_0_0 = get_haplotype_counts(0, het_variants, donor, experiment_key + "1")
    haplotype_counts_0_1 = get_haplotype_counts(0, het_variants, donor, experiment_key + "2")

    
    haplotype_counts_1_0 = get_haplotype_counts(1, het_variants, donor, experiment_key + "1")
    haplotype_counts_1_1 = get_haplotype_counts(1, het_variants, donor, experiment_key + "2")
    
    # Zero counts where evidence of alignment bias
    # Need to expand this to make it less aggressive
    final_counts_0 = []
    final_counts_1 = []
    for a1, a2, a3, b1, b2, b3 in zip(
       haplotype_counts_0, haplotype_counts_0_0, haplotype_counts_0_1,
       haplotype_counts_1, haplotype_counts_1_0, haplotype_counts_1_1
    ):
        if a1 == a2 == a3 and b1 == b2 == b3: # All consistent so nothing to do
            final_counts_0.append(a1)
            final_counts_1.append(b1)
        elif min(a1, a2, a3) >= max(b1, b2, b3): # Take least extreme skew
            final_counts_0.append(min(a1, a2, a3))
            final_counts_1.append(max(b1, b2, b3))
        elif max(a1, a2, a3) <= min(b1, b2, b3): # Take least extreme skew
            final_counts_0.append(max(a1, a2, a3))
            final_counts_1.append(min(b1, b2, b3))
        else: # Data overlaps in unresolvable manner.
            final_counts_0.append(0)
            final_counts_1.append(0)

    return  final_counts_0, final_counts_1


def get_other_haplotype_counts(het_variants, experiment_key):    
    return [v.info.get(experiment_key + "O", 0) for v in het_variants]


def get_het_probabilities_list(phased_variants, donor):
    return [min(get_genotype_probabilities(v, donor)[1], 0.99) for v in phased_variants]


def get_genotype_probabilities(variant, donor):
    genotype_likelihoods = variant.samples[donor].get("PL", None)
    if genotype_likelihoods is None:
        return tuple([0.33] * 3)
    else:
        if len(genotype_likelihoods) != 3:
            print(variant.pos, donor, genotype_likelihoods)
        assert len(genotype_likelihoods) == 3
        
        genotype_likelihoods = list(map(lambda x: 10 ** -x, genotype_likelihoods))
        sum_likelihoods = sum(genotype_likelihoods)
        genotype_likelihoods = [x / sum_likelihoods for x in genotype_likelihoods]
        return tuple(genotype_likelihoods)


def get_phase_probability(var, donor):
    return max(min(1.0 - 10.0 ** (-var.samples[donor].get("PQ", 0.0) / 10.0), 0.99), 0.5)


def get_phase_probabilities(phased_variants, donor):
    return [get_phase_probability(v, donor) for v in phased_variants] 


def get_phasesets(phased_variants, donor):
    return [v.samples[donor].get("PS", None) for v in phased_variants]


def get_variant_accuracy(var):
    # Accuracy from phred is 1 - 10 ** (-P / 10)
    return 1.0 - 10.0 ** (-var.qual / 10.0)


if __name__=="__main__":
    parser = argparse.ArgumentParser(
        prog='Form WASP input files',
        description='Takes donor peak file and count data',
    )
    parser.add_argument("--regions-file", required=True)
    parser.add_argument("--donor-vcf", required=True)
    parser.add_argument("--backbone-vcf", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--experiment-type", required=True)
    parser.add_argument("--has-allele-specific-alignment", default=False, type=bool)
    parser.add_argument("--padding", required=True, type=int)
    parser.add_argument("--max-variant-size", default=1, type=int)
    parser.add_argument("--out")

    args=parser.parse_args()

    backbone_vcf = pysam.VariantFile(args.backbone_vcf)
    donor_vcf = pysam.VariantFile(args.donor_vcf)
    
    # Assumes one donor per VCF
    sample_name = list(donor_vcf.header.samples)[0]
    config = json.load(open(args.config))

    if args.experiment_type not in config:
        print("Missing data for:", args.experiment_type, "for", sample_name)
        sys.exit()
    
    print(config[args.experiment_type])
    experiment_key = config[args.experiment_type]['stub']
    
    stats = pd.read_csv(
        StringIO(pysam.idxstats(config[args.experiment_type]['bam-file'])),
        sep="\t", header=None
    )
    genome_wide_read_count = stats[2].sum()
    print(genome_wide_read_count)
    
    bam_file = pysam.AlignmentFile(config[args.experiment_type]['bam-file'])

    # Files can have different headers. We always have first three columns in this form
    regions = pd.read_csv(args.regions_file, sep="\t", skiprows=0, usecols=range(3))
    regions.columns = ["Chromosome", "Start", "End"]
    output = []

    for _, peak in regions.iterrows():
        # Get the same variants for each donor from the backbone VCF.
        # Most won't exist for this donor. We provide a low genotype score
        # Test variants always come from ATAC peak regions for now
        # This means same test variants for all scenarios
        backbone_variants = get_variants(peak.Chromosome, peak.Start, peak.End, backbone_vcf)

        # Currently skipping multiallelic sites as target SNPs
        # Seriously need to consider examples where this is less complicated
        backbone_variants = [v for v in backbone_variants if len(v.alts) == 1]
        
        # Filter on size of variant
        backbone_variants = [v for v in backbone_variants if max(len(v.ref), len(v.alts[0])) <= args.max_variant_size]
  
        # Compute region_start and region_end
        if peak.End - peak.Start >= 2 * args.padding:
            region_start = peak.Start
            region_end = peak.End
        else:
            mid_point = (peak.Start + peak.End) // 2
            region_start = mid_point - args.padding
            region_end = mid_point + args.padding
        
        # Now extract the evidence variant annotatations.
        # This might allow a padded peak-call for clarity.
        donor_variants = get_variants(peak.Chromosome, region_start, region_end, donor_vcf)
        het_variants = [v for v in donor_variants if set(v.samples[sample_name]["GT"]) == {0, 1}]

        # Currently skipping multiallelic sites as target SNPs
        het_variants = [v for v in het_variants if len(v.alts) == 1]
        
        # Filter on size of variant
        het_variants = [
            v for v in het_variants if max(len(v.ref), len(v.alts[0])) <= args.max_variant_size]
        
        # Compute the read count
        read_count = bam_file.count(peak.Chromosome, region_start, region_end)
        
        # These are the same for all test variants
        hets_positions = [v.pos for v in het_variants]
        hets_genotype_probs = get_het_probabilities_list(het_variants, sample_name)
        hets_phase_probs = get_phase_probabilities(het_variants, sample_name)
        hets_phase_sets = get_phasesets(het_variants, sample_name)
        variant_qualities = list(map(get_variant_accuracy, het_variants))
        no_filters = ["PASS" in v.filter for v in het_variants]

        # Get the haplotype counts from the appropiate alleles
        # Use read counts aligned to the correct genome
        haplotype_counts_0, haplotype_counts_1 = get_haplotypes_counts(
            het_variants, sample_name, experiment_key, args.has_allele_specific_alignment
        )

        other_counts = get_other_haplotype_counts(het_variants, experiment_key)

        for var in backbone_variants:
            # Move this to a function
            genotype = var.samples[sample_name]["GT"]
            # TODO: check that this is the correct way around
            if genotype == (0, 1):
                ref_counts = haplotype_counts_0
                alt_counts = haplotype_counts_1
            elif genotype == (1, 0):
                ref_counts = haplotype_counts_1
                alt_counts = haplotype_counts_0
            elif np.random.randint(2) == 0 and len(set(genotype)) == 1:
                # Follow WASP's approach to randomly asign ref/alt alleles
                # Alternative approach is to provide zeroes. Seems a waste of information
                ref_counts = haplotype_counts_0
                alt_counts = haplotype_counts_1
            elif len(set(genotype)) == 1:
                ref_counts = haplotype_counts_1
                alt_counts = haplotype_counts_0
            else:
                ref_counts = [0 for x in haplotype_counts_1]
                alt_counts = [0 for x in haplotype_counts_0]
                print("Unclear", var.pos, sample_name, genotype)

            matching_vars = [v for v in donor_variants if v.pos == var.pos]
            if len(matching_vars) >= 1:
                if len(matching_vars) > 1:
                    print(
                        "Duplicated variants in file. Take first at pos",
                        [(v.chrom, v.pos, v.ref, v.alts) for v in matching_vars]
                    )

                var_genotype_probs = get_genotype_probabilities(matching_vars[0], sample_name)
                
                genotype_prob = 1.0 * var_genotype_probs[1] + 2.0 * var_genotype_probs[2]
                # TODO: make it so that it matches |,/ from the input
                genotype = "|".join(list(map(str, var.samples[sample_name]["GT"])))
                phase_set = matching_vars[0].samples[sample_name].get("PS", None)
                phase_quality = get_phase_probability(matching_vars[0], sample_name)
                variant_quality = get_variant_accuracy(matching_vars[0])
            else:
                genotype_prob = 0.01
                genotype = "0|0"
                phase_quality = 0.0
                phase_set = None
                variant_quality = 0.01

            genotype = genotype.replace("None", "0")

            # WARNING: WASP uses column indexes so be careful here.
            not_empty = len(het_variants) > 0
            
            adjusted_het_probs = []
            for het_pos, het_prob, linkage_prob, het_phase_set, var_quality, no_filter in zip(
                hets_positions, hets_genotype_probs, hets_phase_probs, hets_phase_sets, variant_qualities, no_filters
            ):
                if het_pos == var.pos:
                    phase_prob = 1.0
                elif phase_set is None or phase_set != het_phase_set:
                    # Don't know phaseset so have 50% chance of being correct
                    phase_prob = 0.5
                else:
                    # Genotype being correct is prob het * linkage * var_prob
#                     phase_prob = linkage_prob * phase_quality
                    phase_prob = 1.0
                
                filter_adjustment = 1.0 if no_filter else 0.5
                adjusted_het_probs.append(het_prob * phase_prob * var_quality * filter_adjustment)

            output.append({
                "CHROM": peak.Chromosome,
                "TEST.SNP.POS": var.pos,
                "TEST.SNP.ID": var.id if var.id else "{}:{}:{}:{}".format(var.chrom, var.pos, var.ref, var.alts[0]),
                "TEST.SNP.REF.ALLELE": var.ref,
                "TEST.SNP.ALT.ALLELE": var.alts[0],
                "TEST.SNP.GENOTYPE": genotype_prob,
                # Ideally this should include possibility not phased for clear understanding
                "TEST.SNP.HAPLOTYPE": genotype,
                "REGION.START": peak.Start,
                "REGION.END": peak.End,
                "REGION.SNP.POS": ";".join(list(map(str, hets_positions))) if not_empty else "NA",
                
                "REGION.SNP.HET.PROB": ";".join(list(map(str, adjusted_het_probs))) if not_empty else "NA",
                "REGION.SNP.LINKAGE.PROB": ";".join(list(map(str, hets_phase_probs))) if not_empty else "NA",
                
                "REGION.SNP.REF.HAP.COUNT": ";".join(list(map(str, ref_counts))) if not_empty else "NA",
                "REGION.SNP.ALT.HAP.COUNT": ";".join(list(map(str, alt_counts))) if not_empty else "NA",
                "REGION.SNP.OTHER.HAP.COUNT": ";".join(list(map(str, other_counts))) if not_empty else "NA",
                
                "REGION.READ.COUNT": read_count,
                "GENOMEWIDE.READ.COUNT": genome_wide_read_count,
            })

    pd.DataFrame(output).to_csv(args.out, sep=" ", index=False, compression='gzip')
