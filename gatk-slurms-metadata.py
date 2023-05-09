# Author: Suhani Balachandran
# Date: Spring 2023
# Project: BSE IW - single cell lm
# Advisors: Mona Singh, Sara Geraghty

# Filename: gatk-slurms-metadata.py
# Purpose: metadata-aware creation of slurm files for running GATK Mutect2.
# Parameters: path to metadata file, path to input bams, path to resources, path to output directory
# Output: slurm files

import pandas as pd
import sys

metadata_path = sys.argv[1]
in_dir = sys.argv[2]
resource_dir = sys.argv[3]
out_dir = sys.argv[4]
slurm_header = "#!/bin/bash\n#SBATCH --ntasks=10 \n#SBATCH--cpus-per-task=2\n#SBATCH --mem 64GB\n#SBATCH --job-name run-gatk\n#SBATCH -o %j.out\n#SBATCH -e %j.err\n#SBATCH --mail-type=END,FAIL\n#SBATCH --mail-user=sb5334@princeton.edu"

all_metadata = pd.read_csv(metadata_path, sep=",")
metadata = all_metadata[["Run", "Cell_type", "neoplastic", "Patient_ID", "Tissue", "tsne_cluster"]]
# did not pre-process properly - removing for now
for x in ['SRR3935146',  'SRR3935082', 'SRR3935528', 'SRR3935616', 'SRR3936026', 'SRR3935755', 'SRR3936282', 'SRR3935995', 'SRR3936407', 'SRR3936504', 'SRR3937007']:
    metadata = metadata[metadata['Run'] != x]

filtered = metadata.sort_values(by=["Patient_ID"])
filtered = filtered[(filtered["neoplastic"] == "Regular") & (filtered["Tissue"] == "Periphery") & (filtered["Cell_type"] == "Oligodendrocyte")]

patient_ids = ["BT_S1", "BT_S2", "BT_S4", "BT_S6"]
patient_refs = [filtered[filtered["Patient_ID"] == x].sample()['Run'].iloc[0] for x in patient_ids]

for p in range(len(patient_ids)):
    count = 0
    subset = metadata[metadata["Patient_ID"] == patient_ids[p]]
    print(subset.head())
    normal_sra = patient_refs[p]

    for x in range(0, len(subset), 8):
        out = open("rungatk-{}-{}.slurm".format(patient_ids[p], count), "w")
        out.write(slurm_header)

        end_index = x+8
        if (x + 8 > len(subset)): end_index = len(subset)

        out.write("\nmodule load GATK")

        for j in range(x, end_index):
            sra = subset['Run'].iloc[j]
            #out.write("\nsrun --ntasks=1 gatk Mutect2 -R {3}hg38.fa -I {0}{1}_final.bam -I {0}{4}_final.bam --germline-resource {3}somatic-hg38_af-only-gnomad.hg38.vcf.gz --panel-of-normals {3}somatic-hg38_1000g_pon.hg38.vcf.gz -O {2}{1}.vcf.gz &".format(in_dir, sra, out_dir, resource_dir, normal_sra))
            out.write("\nsrun --ntasks=1 gatk Mutect2 -R {3}GCF_000001405.40_GRCh38.p14_genomic.fna -I {0}{1}_3.bam -I {0}{4}_3.bam -normal {4} -O {2}{1}.vcf.gz &".format(in_dir, sra, out_dir, resource_dir, normal_sra))

        out.write("\nwait")
        count+=1

    out.close()