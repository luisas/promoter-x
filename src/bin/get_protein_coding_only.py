#!/usr/bin/env python3

import pandas as pd
import sys


gtf = sys.argv[1]
output_gtf = sys.argv[2]
gtf = "/home/luisasantus/Desktop/crg_cluster/projects/PROMOTER-X/src/work/68/c700e0b821f68bac71ad83d66ab913/gencode.v29.annotation.gtf"

gtf_original = pd.read_csv(gtf, sep="\t", comment="#", header=None)
gtf_original["gene_id"] = gtf_original[8].str.extract(r'gene_id "(.*?)";')
gtf_original["transcript_id"] = gtf_original[8].str.extract(r'transcript_id "(.*?)";')
gtf_original = gtf_original[gtf_original[8].str.contains("protein_coding")]
all_protein_coding_genes = gtf_original[gtf_original[2] == "transcript"]
all_protein_coding_genes = all_protein_coding_genes[[0, 3, 4, "transcript_id", 6, "gene_id"]]


# here we want to keep the longes transcript for each gene
# compute the length of each gene
all_protein_coding_genes["length"] = abs(all_protein_coding_genes[4] - all_protein_coding_genes[3])
# for each gene find the longest transcript and keep it only 
all_protein_coding_genes = all_protein_coding_genes.sort_values(by=["gene_id", "length"], ascending=False)
all_protein_coding_genes = all_protein_coding_genes.drop_duplicates(subset="gene_id", keep="first")
#all_protein_coding_genes = all_protein_coding_genes.drop(columns=["length"])

#now just create a list of strings transcript|gene and store in file
all_protein_coding_genes["transcript_gene"] = all_protein_coding_genes["transcript_id"] + "|" + all_protein_coding_genes["gene_id"]
all_protein_coding_genes = all_protein_coding_genes[["transcript_gene"]]
# drop duplicates
all_protein_coding_genes = all_protein_coding_genes.drop_duplicates()
# remove versions from transcript_id and gene_id
all_protein_coding_genes["transcript_gene"] = all_protein_coding_genes["transcript_gene"].str.replace(r'\.\d+', '', regex=True)
all_protein_coding_genes.to_csv(output_gtf, sep="\t", index=False, header  = False)

