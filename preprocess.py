# Author: Suhani Balachandran
# Date: Spring 2023
# Project: BSE IW - single cell lm
# Advisors: Mona Singh, Sara Geraghty

# Filename: preprocess.py
# Purpose: preprocess through GATK best practices workflow
# Parameters: path to list of bams, path to input dir, path to (output) bams dir, path to resource dir
# Output: slurm files

import sys

in_dir = sys.argv[2]
bams_dir = sys.argv[3]
resource_dir = sys.argv[4]
sras = open(sys.argv[1], "r").readlines()

slurm_header = "#!/bin/bash \n#SBATCH --cpus-per-task=2\n#SBATCH --mem 64GB\n#SBATCH --job-name preprocess\n#SBATCH -o %j.out\n#SBATCH -e %j.err\n#SBATCH --mail-type=END,FAIL\n#SBATCH --mail-user=sb5334@princeton.edu\nmodule load GATK"
count = 0
for s in sras:
    sra = s.strip()
    out = open("preprocess-{}.slurm".format(count), "w")
    markdup = "\ngatk MarkDuplicatesSpark -I {0}{1}_final.bam -O {2}{1}_1.bam".format(in_dir, sra, bams_dir)
    fixtags = "\ngatk SetNmMdAndUqTags -R {2}GCF_000001405.40_GRCh38.p14_genomic.fna -I {0}{1}_1.bam -O {0}{1}_2.bam".format(bams_dir, sra, resource_dir)
    baserecal = "\ngatk BaseRecalibrator -I {0}{1}_2.bam -R {2}GCF_000001405.40_GRCh38.p14_genomic.fna --known-sites {2}GCF_000001405.40.gz -O {1}_recal.table".format(bams_dir, sra, resource_dir)
    apply = "\ngatk ApplyBQSR -R {2}GCF_000001405.40_GRCh38.p14_genomic.fna -I {0}{1}_2.bam --bqsr-recal-file {1}_recal.table -O {0}{1}_3.bam".format(bams_dir, sra, resource_dir)

    commands = [slurm_header, markdup, fixtags, baserecal, apply]
    for x in commands: out.write(x)
    out.close()
    count+=1