#!/usr/bin/env python3

import matplotlib.pyplot as plt
import seaborn as sns
import os
import pandas as pd
import seaborn as sns
import glob
import sys
import numpy as np
np.random.seed(42)

promoter_pairs =  pd.read_csv(sys.argv[1])
expression_dir = sys.argv[2]
expression_files = glob.glob(expression_dir+"/**/*.tsv")
outfile = sys.argv[3]




expression_sample_ids = []
expression_dict = {}
for file in expression_files:
    expression = pd.read_csv(file, sep="\t")
    name = file.split(".tsv")[0].split("/")[-1]
    gene_id = expression["gene_id"].astype(str)
    gene_id = gene_id.str.replace(r'\.\d+', '', regex=True)

    complete_name = file.split("/")[-2] + "_" + name
    pme_TPM = expression["pme_TPM"]
    # for each gene_id store the pme_TPM value if it is not already in the dictionary otherwise append the value
    expression_sample_ids.append(complete_name)
    for i in range(len(gene_id)):
        if gene_id[i] in expression_dict:
            expression_dict[gene_id[i]].append(pme_TPM[i])
        else:
            expression_dict[gene_id[i]] = [pme_TPM[i]]


expression = pd.DataFrame(list(expression_dict.items()), columns=['gene_id', 'pme_TPM'])
# add the sample names as columns, it is the same vector for all the genes, in string format
expression["complete_name"] = [list(expression_sample_ids)] * len(expression)

# merge expression data
promoter_pairs_with_expression = promoter_pairs.merge(expression, left_on="gene1", right_on="gene_id")
# keep only the necessary columns
# extract names of first 4 columns
first_columns = promoter_pairs_with_expression.columns[:6]
# append pme_TPM and complete_name
first_columns = [first_columns[0], first_columns[1], first_columns[2], first_columns[3], first_columns[4], first_columns[5], "pme_TPM", "complete_name"]
promoter_pairs_with_expression = promoter_pairs_with_expression[ first_columns ]

second_gene_exp = expression[["gene_id", "pme_TPM"]]
promoter_pairs_with_expression = promoter_pairs_with_expression.merge(second_gene_exp, left_on="gene2", right_on="gene_id")

promoter_pairs_with_expression.reset_index(drop=True, inplace=True)

# save the expression data to a csv file
promoter_pairs_with_expression.to_csv(outfile, index = False)

