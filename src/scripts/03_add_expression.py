import matplotlib.pyplot as plt
import seaborn as sns
import os
import pandas as pd
import seaborn as sns

promoter_pairs = pd.read_csv("./data/promoters/promoter_pairs.csv")
promoter_pairs["Sequence Identity"]
# convert EntrezID to string
promoter_pairs["GeneID1"] = promoter_pairs["GeneID1"].astype(str)
promoter_pairs["GeneID2"] = promoter_pairs["GeneID2"].astype(str)

# get all files in the folder with absolute path
expression_files = os.listdir("./data/expression/liver")
expression_files = [os.path.join("./data/expression/liver", file) for file in expression_files]
# Just read in one expression file and plot the ditribution of the expression values
pc_genes = pd.read_csv("./data/annotation/protein_coding_genes.gtf", sep="\t")
pc_genes_ids = pc_genes["gene_id"]



# READ IN THE EXPERSSIONS, APPEND AND PLOT AND COLOR BY FILE NAME
# read in the expression files

exps = pd.DataFrame()
for file in expression_files:
    expression = pd.read_csv(file, sep="\t")
    name = file.split(".tsv")[0].split("/")[-1]
    print(name)
    expression["name"] = name
    exps = pd.concat([exps, expression], ignore_index=True)

# keep only the protein coding genes
sample_pc = pc_genes_ids.sample(1000)
exps = exps[exps["gene_id"].isin(sample_pc)]
# add pseudocount
exps["TPM"] = exps["TPM"] + 1
# plot histogram and log scale and color by name
sns.histplot(data=exps, x="TPM", hue="name", log_scale=True)
plt.show()






expression = pd.read_csv(expression_files[0], sep="\t")
# get only protein coding genes
expression = expression[expression["gene_id"].isin(pc_genes_ids)]
expression.columns
# extarct only the TPM values and store in a list
exps = expression["TPM"]
# add pseudocount
exps = exps + 1
# plot histogram and log scale
sns.histplot(exps, log_scale=True)
plt.show()



expressions = pd.DataFrame()
for file in expression_files:
    expression = pd.read_csv(file, sep="\t")
    name = file.split(".tsv")[0].split("/")[-1]
    print(name)
    expression["name"] = name
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



# Explore expression range of the genes per file
# plot a histogram of the expression values and color by name
# density plot
sns.histplot(data=expressions, x="pme_TPM_x", hue="name", kde=True)
plt.show()

# check the columns
df_plot = expressions[["GeneID1", "GeneID2", "Sequence Identity", "pme_TPM_x", "pme_TPM_y", "name"]]
df_plot["TPM_similarity"] = abs(df_plot["pme_TPM_x"] - df_plot["pme_TPM_y"])
# make it btw 0 and 1
# df_plot["TPM_similarity"] = df_plot["TPM_similarity"] / df_plot["TPM_similarity"].max()
# check min of TPM_similarity
df_plot["pme_TPM_x"].max()

df_plot
sns.scatterplot(data=df_plot, x="Sequence Identity", y="TPM_similarity", hue= "name")
# put legend outside
plt.show()