# Author: Suhani Balachandran
# Date: Spring 2023
# Project: BSE IW - single cell lm
# Advisors: Mona Singh, Sara Geraghty

# Filename: filtermutect.py
# Purpose: filter variant calls through GATK best practices workflow
# Parameters: path to list of bams, path to input dir, path to (output) bams dir, path to resource dir, path vcf dir
# Output: slurm files

import sys

in_dir = sys.argv[2]
bams_dir = sys.argv[3]
resource_dir = sys.argv[4]
sras = open(sys.argv[1], "r").readlines()

slurm_header = "#!/bin/bash \n#SBATCH --cpus-per-task=1\n#SBATCH --mem 64GB\n#SBATCH --job-name filter\n#SBATCH -o %j.out\n#SBATCH -e %j.err\n#SBATCH --mail-type=END,FAIL\n#SBATCH --mail-user=sb5334@princeton.edu\nmodule load GATK"
count = 0
for s in sras:
    if s not in ['SRR3935146',  'SRR3935082', 'SRR3935528', 'SRR3935616', 'SRR3936026', 'SRR3935755', 'SRR3936282', 'SRR3935995', 'SRR3936407', 'SRR3936504', 'SRR3937007', 'SRR3937390', 'SRR3937566', 'SRR3937857']:
        sra = s.strip()
        out = open("filter-{}.slurm".format(count), "w")

        pileup = "\ngatk GetPileupSummaries -I {0}{1}_3.bam -V {2}GCF_000001405.40.gz -L {2}GCF_000001405.40.gz -O {0}{1}_pileups.table".format(bams_dir, sra, resource_dir)
        contamination = "\ngatk CalculateContamination -I {0}{1}_pileups.table -O {0}{1}_contamination.table".format(bams_dir, sra)
        collect = "\ngatk CollectF1R2Counts -R {2}GCF_000001405.40_GRCh38.p14_genomic.fna -I {0}{1}_3.bam -O {0}{1}_prior.tsv --tumor-segmentation {0}{1}_segments.tsv".format(bams_dir, sra, resource_dir)
        tar = "\ntar -czvf {0}{1}_f1r2.tar.gz {0}{1}_prior.tsv".format(bams_dir, sra)
        learn = "\ngatk LearnReadOrientationModel {0}{1}_f1r2.tar.gz -O {0}{1}_artifact-prior.tar.gz".format(bams_dir, sra)
        filtered = "\ngatk FilterMutectCalls -R {2}GCF_000001405.40_GRCh38.p14_genomic.fna -V {3}{1}.vcf.gz --contamination-table {0}{1}_contamination.table --tumor-segmentation {0}{1}_segments.tsv -O {3}{1}_filtered.vcf.gz".format(bams_dir, sra, resource_dir, sys.argv[5])

        commands = [slurm_header, pileup, contamination, collect, tar, learn, filtered]
        for x in commands: out.write(x)
        out.close()
        count+=1
