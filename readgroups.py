# Author: Suhani Balachandran
# Date: Spring 2023
# Project: BSE IW - single cell lm
# Advisors: Mona Singh, Sara Geraghty

# Filename: readgroups.py
# Purpose: add read groups to files before calling GATK
# Parameters: path to metadata file, path to input bams, path to output
# Output: slurm files

import pandas as pd
import sys

# GATK TERMS --> MY TERMS
# FLOWCELL_BARCODE --> PLATE.ID
# LANE --> WELL
# SAMPLE_BARCODE --> SAMPLE_NAME
# SAMPLE --> RUN
# LIBRARY IDENTIFIER --> EXPERIMENT

metadata_path = sys.argv[1]
in_dir = sys.argv[2]
out_dir = sys.argv[3]
slurm_header = "#!/bin/bash\n#SBATCH --ntasks=10 \n#SBATCH--cpus-per-task=2\n#SBATCH --mem 64GB\n#SBATCH --job-name add-rg\n#SBATCH -o %j.out\n#SBATCH -e %j.err\n#SBATCH --mail-type=END,FAIL\n#SBATCH --mail-user=sb5334@princeton.edu"

metadata = pd.read_csv(metadata_path, sep=",")

count = 0
for x in range(0, len(metadata), 10):
    out = open("addrg-{}.slurm".format(count), "w")
    out.write(slurm_header)
    out.write("\nmodule load GATK")

    end_index = x + 10
    if (x + 10 > len(metadata)): end_index = len(metadata)

    for y in range(x, end_index):
        row = metadata.iloc[y]
        id = "{}.{}".format(row['Plate.ID'], row['Well'])
        lb = "{}".format(row['Experiment'])
        pl = "{}".format(row['Platform'])
        pu = "{}.{}.{}".format(row['Plate.ID'], row['Well'], row['Sample Name'])
        sm = "{}".format(row['Run'])
        
        out.write("\nsrun --ntasks=1 gatk AddOrReplaceReadGroups -I {0}{1}_Aligned.out.bam -O {2}{1}_final.bam -RGID {3} -RGLB {4} -RGPL {5} -RGPU {6} -RGSM {7} &".format(in_dir, row['Run'], out_dir, id, lb, pl, pu, sm))
    
    out.write("\nwait")
    out.close()
    count+=1
