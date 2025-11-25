import json
import os
import argparse


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='SkewConfig',
        description='Takes Donor + bam file patterns. Builds config for analysis',
    )
    parser.add_argument("--donor")
    parser.add_argument("--out")

    args=parser.parse_args()

    new_human_data = "/project/Wellcome_Discovery/shared/01_in_house_data/01_50Donors"
    
    # RNA-seq ref alignment
    RNA_SEQ_PATTERN = os.path.join(new_human_data,                                 "03_RNASeq/01_CATCH_ALL/Donor50_ref/results/07_indices/{donor}_PolyAplus_Ery_d13_CS_hg38.bam"
    )
    RNA_SEQ_PATTERN_FORWARD = os.path.join(new_human_data,                                           "03_RNASeq/02_STRANDED/Donor50_ref/{donor}_PolyAplus_Ery_d13_CS_hg38_forward.bam"
    )
    RNA_SEQ_PATTERN_REVERSE = os.path.join(new_human_data,                                           "03_RNASeq/02_STRANDED/Donor50_ref/{donor}_PolyAplus_Ery_d13_CS_hg38_reverse.bam"
    )

    # RNA-seq allele1 alignment
    RNA_SEQ_1_PATTERN = os.path.join(new_human_data,                                   "03_RNASeq/01_CATCH_ALL/Donor50_allele1/results/07_indices/{donor}_PolyAplus_Ery_d13_CS_hg38.bam"
    )
    RNA_SEQ_1_PATTERN_FORWARD = os.path.join(new_human_data,                                           "03_RNASeq/02_STRANDED/Donor50_allele1/{donor}_PolyAplus_Ery_d13_CS_hg38_forward.bam"
    )
    RNA_SEQ_1_PATTERN_REVERSE = os.path.join(new_human_data,                                           "03_RNASeq/02_STRANDED/Donor50_allele1/{donor}_PolyAplus_Ery_d13_CS_hg38_reverse.bam"
    )

    # RNA-seq allele2 alignment
    RNA_SEQ_2_PATTERN = os.path.join(new_human_data,                                   "03_RNASeq/01_CATCH_ALL/Donor50_allele2/results/07_indices/{donor}_PolyAplus_Ery_d13_CS_hg38.bam"
    )
    RNA_SEQ_2_PATTERN_FORWARD = os.path.join(new_human_data,                                           "03_RNASeq/02_STRANDED/Donor50_allele2/{donor}_PolyAplus_Ery_d13_CS_hg38_forward.bam"
    )
    RNA_SEQ_2_PATTERN_REVERSE = os.path.join(new_human_data,                                           "03_RNASeq/02_STRANDED/Donor50_allele2/{donor}_PolyAplus_Ery_d13_CS_hg38_reverse.bam"
    )

    # PMINUS ref alignment
    RNA_SEQ_PMINUS_PATTERN = os.path.join(new_human_data,                                   "03_RNASeq/01_CATCH_ALL/Donor50_ref/results/07_indices/{donor}_Nascent_Ery_d13_CS_hg38.bam"
    )
    RNA_SEQ_PMINUS_PATTERN_FORWARD = os.path.join(new_human_data,                                 "03_RNASeq/02_STRANDED/Donor50_ref/{donor}_Nascent_Ery_d13_CS_hg38_forward.bam"
    )
    RNA_SEQ_PMINUS_PATTERN_REVERSE = os.path.join(new_human_data,                                       "03_RNASeq/02_STRANDED/Donor50_ref/{donor}_Nascent_Ery_d13_CS_hg38_reverse.bam"
    )
    
    # PMINUS allele1 alignment
    RNA_SEQ_PMINUS_1_PATTERN = os.path.join(new_human_data,                                   "03_RNASeq/01_CATCH_ALL/Donor50_allele1/results/07_indices/{donor}_Nascent_Ery_d13_CS_hg38.bam"
    )
    RNA_SEQ_PMINUS_1_PATTERN_FORWARD = os.path.join(new_human_data,                                       "03_RNASeq/02_STRANDED/Donor50_allele1/{donor}_Nascent_Ery_d13_CS_hg38_forward.bam"
    )
    RNA_SEQ_PMINUS_1_PATTERN_REVERSE = os.path.join(new_human_data,                                       "03_RNASeq/02_STRANDED/Donor50_allele1/{donor}_Nascent_Ery_d13_CS_hg38_reverse.bam"
    )

    # PMINUS allele2 alignment
    RNA_SEQ_PMINUS_2_PATTERN = os.path.join(new_human_data,                                   "03_RNASeq/01_CATCH_ALL/Donor50_allele2/results/07_indices/{donor}_Nascent_Ery_d13_CS_hg38.bam"
    )
    RNA_SEQ_PMINUS_2_PATTERN_FORWARD = os.path.join(new_human_data,                                       "03_RNASeq/02_STRANDED/Donor50_allele2/{donor}_Nascent_Ery_d13_CS_hg38_forward.bam"
    )
    RNA_SEQ_PMINUS_2_PATTERN_REVERSE = os.path.join(new_human_data,                                       "03_RNASeq/02_STRANDED/Donor50_allele2/{donor}_Nascent_Ery_d13_CS_hg38_reverse.bam"
    )
    
    ATAC_PATTERN = os.path.join(new_human_data, "01_ATAC/01_CATCH_ALL/Donor50_ref/results/07_merge/{donor}_hg38.bam")
    ATAC_1_PATTERN = os.path.join(new_human_data, "01_ATAC/01_CATCH_ALL/Donor50_allele1/results/07_merge/{donor}_hg38.bam")
    ATAC_2_PATTERN = os.path.join(new_human_data, "01_ATAC/01_CATCH_ALL/Donor50_allele2/results/07_merge/{donor}_hg38.bam") 
    
    CTCF_PATTERN = os.path.join(new_human_data, "02_CTCF/01_CATCH_ALL/Donor50_CTCF_ref/results/07_merge/{donor}_hg38.bam")
    CTCF_1_PATTERN = os.path.join(new_human_data, "02_CTCF/01_CATCH_ALL/Donor50_CTCF_allele1/results/07_merge/{donor}_hg38.bam")
    CTCF_2_PATTERN = os.path.join(new_human_data, "02_CTCF/01_CATCH_ALL/Donor50_CTCF_allele2/results/07_merge/{donor}_hg38.bam")

    RAD21_PATTERN = os.path.join(new_human_data, "04_RAD21/01_CATCH_ALL/RAD21_ref/results/07_merge/{donor}_hg38.bam")
    RAD21_1_PATTERN = os.path.join(new_human_data, "04_RAD21/01_CATCH_ALL/RAD21_allele1/results/07_merge/{donor}_hg38.bam")
    RAD21_2_PATTERN = os.path.join(new_human_data, "04_RAD21/01_CATCH_ALL/RAD21_allele2/results/07_merge/{donor}_hg38.bam")
    
    H3K27ac_PATTERN = os.path.join(new_human_data, "05_H3K27ac/01_CATCH_ALL/Donor50_H3K27ac_ref/results/07_merge/{donor}_hg38.bam"
    )
    H3K27ac_1_PATTERN = os.path.join(new_human_data, "05_H3K27ac/01_CATCH_ALL/Donor50_H3K27ac_allele1/results/07_merge/{donor}_hg38.bam"
    )
    H3K27ac_2_PATTERN = os.path.join(new_human_data, "05_H3K27ac/01_CATCH_ALL/Donor50_H3K27ac_allele2/results/07_merge/{donor}_hg38.bam"
    )

    H3K4me1_PATTERN = os.path.join(new_human_data, "06_H3K4me1/01_CATCH_ALL/Donor50_H3K4me1_ref/results/07_merge/{donor}_hg38.bam"
    )
    H3K4me1_1_PATTERN = os.path.join(new_human_data, "06_H3K4me1/01_CATCH_ALL/Donor50_H3K4me1_allele1/results/07_merge/{donor}_hg38.bam"
    )
    H3K4me1_2_PATTERN = os.path.join(new_human_data, "06_H3K4me1/01_CATCH_ALL/Donor50_H3K4me1_allele2/results/07_merge/{donor}_hg38.bam"
    )

    H3K4me3_PATTERN = os.path.join(new_human_data, "07_H3K4me3/01_CATCH_ALL/Donor50_H3K4me3_ref/results/07_merge/{donor}_hg38.bam"
    )
    H3K4me3_1_PATTERN = os.path.join(new_human_data, "07_H3K4me3/01_CATCH_ALL/Donor50_H3K4me3_allele1/results/07_merge/{donor}_hg38.bam"
    )
    H3K4me3_2_PATTERN = os.path.join(new_human_data, "07_H3K4me3/01_CATCH_ALL/Donor50_H3K4me3_allele2/results/07_merge/{donor}_hg38.bam"
    )

    
    liftover_allele1 = os.path.join(new_human_data, "00_WGS/00_10x_GATK/02_phased_genomes/{donor}.1.chain")
    liftover_allele2 = os.path.join(new_human_data, "00_WGS/00_10x_GATK/02_phased_genomes/{donor}.2.chain")
    
    config={
        "ATAC": {
            "bam-file": ATAC_PATTERN.format(donor=args.donor),
            "stub": "A"
        },
        "ATAC_allele1": {
            "bam-file": ATAC_1_PATTERN.format(donor=args.donor),
            "liftover-file": liftover_allele1.format(donor=args.donor),
            "stub": "A1"
        },
        "ATAC_allele2": {
            "bam-file": ATAC_2_PATTERN.format(donor=args.donor),
            "liftover-file": liftover_allele2.format(donor=args.donor),
            "stub": "A2"
        },
        "CTCF": {
            "bam-file": CTCF_PATTERN.format(donor=args.donor),
            "stub": "C"
        },
        "CTCF_allele1": {
            "bam-file": CTCF_1_PATTERN.format(donor=args.donor),
            "liftover-file": liftover_allele1.format(donor=args.donor),
            "stub": "C1"
        },
        "CTCF_allele2": {
            "bam-file": CTCF_2_PATTERN.format(donor=args.donor),
            "liftover-file": liftover_allele2.format(donor=args.donor),
            "stub": "C2"
        },
        # All RNA-seq on merged
        "RNASeq": {
            "bam-file": RNA_SEQ_PATTERN.format(donor=args.donor),
            "stub": "R"
        },
        "RNASeq_allele1": {
            "bam-file": RNA_SEQ_1_PATTERN.format(donor=args.donor),
            "liftover-file": liftover_allele1.format(donor=args.donor),
            "stub": "R1"
        },
        "RNASeq_allele2": {
            "bam-file": RNA_SEQ_2_PATTERN.format(donor=args.donor),
            "liftover-file": liftover_allele2.format(donor=args.donor),
            "stub": "R2"
        },
        "PolyMinus": {
            "bam-file": RNA_SEQ_PMINUS_PATTERN.format(donor=args.donor),
            "stub": "M"
        },
        "PolyMinus_allele1": {
            "bam-file": RNA_SEQ_PMINUS_1_PATTERN.format(donor=args.donor),
            "liftover-file": liftover_allele1.format(donor=args.donor),
            "stub": "M1"
        },
        "PolyMinus_allele2": {
            "bam-file": RNA_SEQ_PMINUS_2_PATTERN.format(donor=args.donor),
            "liftover-file": liftover_allele2.format(donor=args.donor),
            "stub": "M2"
        },
        # All RNA-seq on forward strand
        "RNASeq_forward": {
            "bam-file": RNA_SEQ_PATTERN_FORWARD.format(donor=args.donor),
            "stub": "RF"
        },
        "RNASeq_allele1_forward": {
            "bam-file": RNA_SEQ_1_PATTERN_FORWARD.format(donor=args.donor),
            "liftover-file": liftover_allele1.format(donor=args.donor),
            "stub": "RF1"
        },
        "RNASeq_allele2_forward": {
            "bam-file": RNA_SEQ_2_PATTERN_FORWARD.format(donor=args.donor),
            "liftover-file": liftover_allele2.format(donor=args.donor),
            "stub": "RF2"
        },
        "PolyMinus_forward": {
            "bam-file": RNA_SEQ_PMINUS_PATTERN_FORWARD.format(donor=args.donor),
            "stub": "MF"
        },
        "PolyMinus_allele1_forward": {
            "bam-file": RNA_SEQ_PMINUS_1_PATTERN_FORWARD.format(donor=args.donor),
            "liftover-file": liftover_allele1.format(donor=args.donor),
            "stub": "MF1"
        },
        "PolyMinus_allele2_forward": {
            "bam-file": RNA_SEQ_PMINUS_2_PATTERN_FORWARD.format(donor=args.donor),
            "liftover-file": liftover_allele2.format(donor=args.donor),
            "stub": "MF2"
        },
        # All RNA-seq on reverse strand
        "RNASeq_reverse": {
            "bam-file": RNA_SEQ_PATTERN_REVERSE.format(donor=args.donor),
            "stub": "RR"
        },
        "RNASeq_allele1_reverse": {
            "bam-file": RNA_SEQ_1_PATTERN_REVERSE.format(donor=args.donor),
            "liftover-file": liftover_allele1.format(donor=args.donor),
            "stub": "RR1"
        },
        "RNASeq_allele2_reverse": {
            "bam-file": RNA_SEQ_2_PATTERN_REVERSE.format(donor=args.donor),
            "liftover-file": liftover_allele2.format(donor=args.donor),
            "stub": "RR2"
        },
        "PolyMinus_reverse": {
            "bam-file": RNA_SEQ_PMINUS_PATTERN_REVERSE.format(donor=args.donor),
            "stub": "MR"
        },
        "PolyMinus_allele1_reverse": {
            "bam-file": RNA_SEQ_PMINUS_1_PATTERN_REVERSE.format(donor=args.donor),
            "liftover-file": liftover_allele1.format(donor=args.donor),
            "stub": "MR1"
        },
        "PolyMinus_allele2_reverse": {
            "bam-file": RNA_SEQ_PMINUS_2_PATTERN_REVERSE.format(donor=args.donor),
            "liftover-file": liftover_allele2.format(donor=args.donor),
            "stub": "MR2"
        },
        
        "H3K27ac": {
            "bam-file": H3K27ac_PATTERN.format(donor=args.donor),
            "stub": "K27"
        },
        "H3K27ac_allele1": {
            "bam-file": H3K27ac_1_PATTERN.format(donor=args.donor),
            "liftover-file": liftover_allele1.format(donor=args.donor),
            "stub": "K271"
        },
        "H3K27ac_allele2": {
            "bam-file": H3K27ac_2_PATTERN.format(donor=args.donor),
            "liftover-file": liftover_allele2.format(donor=args.donor),
            "stub": "K272"
        },
        "H3K4me1": {
            "bam-file": H3K4me1_PATTERN.format(donor=args.donor),
            "stub": "K41"
        },
        "H3K4me1_allele1": {
            "bam-file": H3K4me1_1_PATTERN.format(donor=args.donor),
            "liftover-file": liftover_allele1.format(donor=args.donor),
            "stub": "K411"
        },
        "H3K4me1_allele2": {
            "bam-file": H3K4me1_2_PATTERN.format(donor=args.donor),
            "liftover-file": liftover_allele2.format(donor=args.donor),
            "stub": "K412"
        },
        "H3K4me3": {
            "bam-file": H3K4me3_PATTERN.format(donor=args.donor),
            "stub": "K43"
        },
        "H3K4me3_allele1": {
            "bam-file": H3K4me3_1_PATTERN.format(donor=args.donor),
            "liftover-file": liftover_allele1.format(donor=args.donor),
            "stub": "K431"
        },
        "H3K4me3_allele2": {
            "bam-file": H3K4me3_2_PATTERN.format(donor=args.donor),
            "liftover-file": liftover_allele2.format(donor=args.donor),
            "stub": "K432"
        },
        "RAD21": {
            "bam-file": RAD21_PATTERN.format(donor=args.donor),
            "stub": "D",
        },
        "RAD21_allele1": {
            "bam-file": RAD21_1_PATTERN.format(donor=args.donor),
            "liftover-file": liftover_allele1.format(donor=args.donor),
            "stub": "D1",
        },
        "RAD21_allele2": {
            "bam-file": RAD21_2_PATTERN.format(donor=args.donor),
            "liftover-file": liftover_allele2.format(donor=args.donor),
            "stub": "D2",
        },
    }
    # Add d7, and d10 data where appropriate
    old_config = config.copy()
    for key, value in old_config.items():
        if "d13" in value["bam-file"]:
            for day in ["d7", "d10"]:
                new_key = key + "_" + day
                
                new_value = value.copy()
                new_value["bam-file"] = new_value["bam-file"].replace("d13", day)
                new_value["stub"] = day + value["stub"]
                
                config[new_key] = new_value
    
    dict_keys=list(config.keys())
    for key in dict_keys:
        if not os.path.exists(config[key]["bam-file"]):
            print("Missing {} data for {}".format(key, args.donor))
            print("Expected file: {}".format(repr(config[key]["bam-file"])))
            del config[key] 
    
    with open(args.out, "w") as fp:
        json.dump(config, fp, indent=4)
