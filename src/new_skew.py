import pysam
import argparse
from collections import Counter
import json
from pyliftover import LiftOver
import gc


def get_alleles(bam_file, liftover_file, chrom, pos, padding=2, min_base_quality=10, min_mapping_quality=10):
    def padding_filter(r, adjusted_position):
        # Don't trust when position too close to end of read
        return padding <= adjusted_position - 1 - r.reference_start < len(r.query_sequence) - padding

    alleles = []
    indels = []
    if liftover_file is not None:
        liftovers = liftover_file.convert_coordinate(chrom, pos, "-")
        if len(liftovers) != 1:
            # What to do here? How common is it?
            # Ideally express that the liftover failed
            return dict(Counter(alleles)), dict(Counter(indels))
        
        adjusted_position = liftovers[0][1]
    else:
        adjusted_position = pos
    
    for pileupcolumn in bam_file.pileup(
        chrom,
        adjusted_position - 1,
        adjusted_position,
        min_base_quality=min_base_quality,
        min_mapping_quality=min_mapping_quality,
    ):
        if pileupcolumn.pos == adjusted_position - 1:
            for pileupread in pileupcolumn.pileups:
                if not padding_filter(pileupread.alignment, adjusted_position):
                    continue
                
                indels.append(pileupread.indel)
                
                if not pileupread.is_del and not pileupread.is_refskip:
                    alleles.append(pileupread.alignment.query_sequence[pileupread.query_position])
                else:
                    alleles.append(None)

    return dict(Counter(alleles)), dict(Counter(indels))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
                prog='NewSkewPipeline',
                description='Takes input VCF. Finds counts at sites in VCF from bam(s) files',
    )
    parser.add_argument("--input-vcf", required=True)

    # TODO: Extend for 
    parser.add_argument("--chrom", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--key", default=None)
    parser.add_argument("--out", required=True)

    args=parser.parse_args()

    variant_file = pysam.VariantFile(args.input_vcf)
    
    with open(args.config) as fp:
        bam_file_config = json.load(fp)
    
    # Only process one bam at a time
    if args.key is not None:
        bam_file_config = {
            args.key: bam_file_config[args.key]
        }
    
    # Set the files
    for bam_key in bam_file_config.keys():
        bam_file_config[bam_key]["AlignmentFile"] = pysam.AlignmentFile(
            bam_file_config[bam_key]["bam-file"]
        )
        if "liftover-file" in bam_file_config[bam_key]:
            bam_file_config[bam_key]["LiftOver"] = LiftOver(bam_file_config[bam_key]["liftover-file"])
            # We only need LiftOver for one contig. Clean up
            del bam_file_config[bam_key]["LiftOver"].chain_file.chains
            bam_file_config[bam_key]["LiftOver"].chain_file.chain_index = {
                args.chrom: bam_file_config[bam_key]["LiftOver"].chain_file.chain_index[args.chrom]
            }
            gc.collect()
        else:
            bam_file_config[bam_key]["LiftOver"] = None

    new_header = pysam.VariantHeader()
    for line in str(variant_file.header).split("\n"):
        if "##" not in line:
            continue
        if "##INFO" in line:
            continue
        if "##FORMAT" in line:
            continue
        else:
            new_header.add_line(line)
    
    for bam_key in bam_file_config.keys():
        stub = bam_file_config[bam_key]["stub"]
        new_header.add_line(
            "##INFO=<ID={}R,Number=1,Type=Integer,Description=\"Reference allele observed count for {}\">".format(
                stub, bam_key)
        )
        new_header.add_line(
            "##INFO=<ID={}A,Number=1,Type=Integer,Description=\"Alternative allele observed count for {}\">".format(
                stub, bam_key)
        )
        new_header.add_line(
            "##INFO=<ID={}O,Number=1,Type=Integer,Description=\"Other allele observed count for {}\">".format(
                stub, bam_key)
        )
    
    new_header.add_line(
        "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">"
    )
    new_header.add_line("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")
    new_header.add_line(
        "##FORMAT=<ID=GQ,Number=1,Type=Float,Description=\"Genotype Quality\">"
    )
    new_header.add_line(
        "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification\">"
    )
    new_header.add_line("##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"ID of Phase Set for Variant\">")
    new_header.add_line(
        "##FORMAT=<ID=PQ,Number=1,Type=Integer,Description=\"Phred QV indicating probability at this variant is incorrectly phased\">"
    )
    new_header.add_samples(variant_file.header.samples)

    samples = list(new_header.samples)

    with pysam.VariantFile(args.out, "w", header=new_header) as output_vcf:
        for variant in variant_file.fetch(contig=args.chrom):
            info=None
            if len(variant.alts) == 1:
                indel_length = len(variant.alts[0]) - len(variant.ref)

                info = {}
                for bam_key in bam_file_config.keys():
                    # TODO: Run with different quality filters
                    alleles, indels = get_alleles(
                        bam_file_config[bam_key]["AlignmentFile"],
                        bam_file_config[bam_key]["LiftOver"],
                        variant.chrom,
                        variant.pos,
                        padding=2, min_base_quality=10, min_mapping_quality=10
                    )
                    
                    stub = bam_file_config[bam_key]["stub"]
                    
                    if indel_length == 0:
                        # Count based upon alleles at site.
                        info[stub + "R"] = alleles.get(variant.ref, 0)
                        info[stub + "A"] = alleles.get(variant.alts[0], 0)
                        total_ref_alt = info[stub + "R"] + info[stub + "A"]
                        info[stub + "O"] = sum(alleles.values()) - total_ref_alt
                    else:
                        # Count based upon indels after site.
                        info[stub + "R"] = indels.get(0, 0)
                        info[stub + "A"] = indels.get(indel_length, 0)
                        total_ref_alt = info[stub + "R"] + info[stub + "A"]
                        info[stub + "O"] = sum(indels.values()) - total_ref_alt
                        
            new_record = output_vcf.new_record(
                contig = variant.contig, start=variant.start, stop=variant.stop, alleles=variant.alleles,
                id=variant.id, qual=variant.qual, filter=variant.filter,
                info=info
            )
            for sample in samples:
                new_record.samples[sample]["GT"] = variant.samples[sample]["GT"]
                if variant.samples[sample].phased:
                    new_record.samples[sample].phased = True
                for key in ["GQ", "PL", "AD", "PS", "PQ"]:
                    if key in variant.samples[sample]:
                        new_record.samples[sample][key] = variant.samples[sample][key]
            
            output_vcf.write(new_record)
