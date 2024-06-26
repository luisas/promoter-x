#!/usr/bin/env python3
import sys
import pandas as pd

sw_file = sys.argv[1]
promoter_lenght = int(sys.argv[2])
outfile = sys.argv[3]


def sw_to_sim(sw_file, promoter_lenght):
    df = pd.read_csv(sw_file, sep='\t', header = None, on_bad_lines='skip')
    df.columns = ["promoter_similarity_sw", "gene1", "gene2"]
    # remove any spaces in the gene names
    df["gene1"] = df["gene1"].str.strip()
    df["gene2"] = df["gene2"].str.strip()
    df["promoter_similarity_sw"] =  (df["promoter_similarity_sw"]*100)/promoter_lenght
    return df

scores = sw_to_sim(sw_file, promoter_lenght)

scores.to_csv(outfile, sep=',', index=False, header=True)