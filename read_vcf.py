# Author: Suhani Balachandran
# Date: Spring 2023
# Project: BSE IW - single cell lm
# Advisors: Mona Singh, Sara Geraghty

# Filename: read_vcf.py
# Purpose: process VCFs to determine if there are mutations in the driver gene areas
# Parameters: path to driver locations file, path to metadata file, path to input vcfs, path to output CSV
# Output: mutation status CSV (sample x driver gene)

import io
import sys
import os
import pandas as pd

driver_locs = sys.argv[1]
metadata_path = sys.argv[2]
input_dir = sys.argv[3]
output_dir = sys.argv[4]

# FUNCTION FROM https://gist.github.com/dceoy/99d976a2c01e7f0ba1c813778f9db744
def read_vcf(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})

all_vcfs = [x for x in os.listdir(input_dir) if x.endswith('vcf')]
driver_df = pd.read_csv(driver_locs, sep=",")
out_df = pd.DataFrame(columns= (['Sample'] + driver_df['GENE']))
metadata = pd.read_csv(metadata_path, sep=",")

for vcf in all_vcfs:
    df = read_vcf("{}{}".format(input_dir, vcf))
    meta = metadata[metadata['Run'] == vcf.split("_")[0]]

    row = {'Sample': "{}.{}".format(meta['Plate.ID'], meta['Well'])} # match id name in provided count matrix
    for d in len(driver_df):
        mut = df[df['CHROM'] == driver_df['CHROM'].iloc[d] & df['POS'] >= driver_df['START'].iloc[d] & df['POS'] <= driver_df['END'].iloc[d]]
        if len(mut) > 0: row[driver_df['GENE'].iloc[d]] = 1
        else: row[driver_df['GENE'].iloc[d]] = 0
    out_df.append(pd.DataFrame(row), ignore_index=True)

result = open("mutations.csv", "w")
out_df.to_csv(result, index=False)
result.close()