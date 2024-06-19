import matplotlib.pyplot as plt
import seaborn as sns
import os
import pandas as pd
import seaborn as sns
import glob
import numpy as np

promoter_pairs = pd.read_csv("./data/promoters/promoter_pairs_test.csv")
promoter_pairs["Sequence Identity"]
# plot histogram of the sequence identity
sns.histplot(promoter_pairs["Sequence Identity"])
plt.show()

# count how many promoters we have as GeneID1 and GeneID2
promoter_pairs["GeneID1"].value_counts()
promoter_pairs["GeneID2"].value_counts()


# convert EntrezID to string
promoter_pairs["GeneID1"] = promoter_pairs["GeneID1"].astype(str)
promoter_pairs["GeneID2"] = promoter_pairs["GeneID2"].astype(str)

# get all files in the folder with absolute path
expression_files = glob.glob("./data/expression/*/*.tsv")
# Just read in one expression file and plot the ditribution of the expression values
pc_genes = pd.read_csv("./data/annotation/protein_coding_genes.gtf", sep="\t")
pc_genes_ids = pc_genes["gene_id"]



expressions = pd.DataFrame()
for file in expression_files:
    expression = pd.read_csv(file, sep="\t")
    name = file.split(".tsv")[0].split("/")[-1]
    print(name)
    expression["name"] = name
    expression["tissue"] = file.split("/")[-2]
    expression["complete_name"] = file.split("/")[-2] + "_" + name
    # load promoter pairs
    expression["gene_id"] = expression["gene_id"].astype(str)

    # merge expression data
    promoter_pairs_with_expression = promoter_pairs.merge(expression, left_on="GeneID1", right_on="gene_id")
    # keep only the necessary columns
    promoter_pairs_with_expression = promoter_pairs_with_expression[[ "name", "GeneID1", "GeneID2", "Sequence Identity", "pme_TPM"]]


    second_gene_exp = expression[["gene_id", "pme_TPM"]]
    promoter_pairs_with_expression = promoter_pairs_with_expression.merge(second_gene_exp, left_on="GeneID2", right_on="gene_id")

    promoter_pairs_with_expression.reset_index(drop=True, inplace=True)
    expressions = pd.concat([expressions, promoter_pairs_with_expression], ignore_index=True)



expressions
# Explore expression range of the genes per file
# plot a histogram of the expression values and color by name
# density plot
sns.histplot(data=expressions, x="pme_TPM_y", hue="name", kde=True, log_scale=True)
plt.show()

# check the columns
df_plot = expressions[["GeneID1", "GeneID2", "Sequence Identity", "pme_TPM_x", "pme_TPM_y", "name"]]
df_plot["TPM_similarity"] = abs(df_plot["pme_TPM_x"] - df_plot["pme_TPM_y"])



# log the TPM_similarity
df_plot["TPM_similarity"] = df_plot["TPM_similarity"] + 1
df_plot["TPM_similarity"] = df_plot["TPM_similarity"].apply(lambda x: np.log(x))
sns.scatterplot(data=df_plot, x="Sequence Identity", y="TPM_similarity", hue= "complete_name")
# y label 
plt.ylabel("log(TPM similarity)")
# put legend outside
plt.show()