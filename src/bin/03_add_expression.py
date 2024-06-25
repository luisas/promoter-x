import matplotlib.pyplot as plt
import seaborn as sns
import os
import pandas as pd
import seaborn as sns
import glob
import numpy as np
np.random.seed(42)
blue_palette = sns.color_palette("Blues", n_colors=10)
green_palette = sns.color_palette("Reds", n_colors=10)
# Create a custom palette with the first two shades of blue and the second two shades of green
custom_palette = [blue_palette[6], blue_palette[9], green_palette[6], green_palette[9]]
# Set the palette
sns.set_palette(custom_palette)



promoter_pairs = pd.read_csv("./data/promoters/promoter_pairs_test_FULL.csv")
expression_files = glob.glob("./data/expression/*/*.tsv")
pc_genes = pd.read_csv("./data/annotation/protein_coding_genes.gtf", sep="\t")
pc_genes_ids = pc_genes["gene_id"]



expressions = pd.DataFrame()
for file in expression_files:
    expression = pd.read_csv(file, sep="\t")
    name = file.split(".tsv")[0].split("/")[-1]
    expression["name"] = name
    expression["tissue"] = file.split("/")[-2]
    expression["complete_name"] = file.split("/")[-2] + "_" + name
    # load promoter pairs
    expression["gene_id"] = expression["gene_id"].astype(str)

    # merge expression data
    promoter_pairs_with_expression = promoter_pairs.merge(expression, left_on="GeneID1", right_on="gene_id")
    # keep only the necessary columns
    promoter_pairs_with_expression = promoter_pairs_with_expression[[ "name", "GeneID1", "GeneID2", "Sequence Identity", "pme_TPM", "complete_name"]]


    second_gene_exp = expression[["gene_id", "pme_TPM"]]
    promoter_pairs_with_expression = promoter_pairs_with_expression.merge(second_gene_exp, left_on="GeneID2", right_on="gene_id")

    promoter_pairs_with_expression.reset_index(drop=True, inplace=True)
    expressions = pd.concat([expressions, promoter_pairs_with_expression], ignore_index=True)


# check the columns
df_plot = expressions[["GeneID1", "GeneID2", "Sequence Identity", "pme_TPM_x", "pme_TPM_y", "complete_name"]]
df_plot = df_plot[df_plot["complete_name"] == "liver_ENCFF658MFL"]
df_plot = df_plot[df_plot["GeneID1"] != df_plot["GeneID2"]]
df_plot["pme_TPM_x"] = df_plot["pme_TPM_x"] + 0.01
df_plot["pme_TPM_y"] = df_plot["pme_TPM_y"] + 0.01

promoters_seq_id_1 = df_plot[df_plot["Sequence Identity"] == 1]
promoters_seq_id_1
# plot the pme_TPM_x vs pme_TPM_y
sns.scatterplot(data=promoters_seq_id_1, x="pme_TPM_x", y="pme_TPM_y", hue="Sequence Identity", s = 3)
plt.show()

# now 



df_plot_nozeros = df_plot[df_plot["pme_TPM_x"] > 0]
df_plot_nozeros = df_plot[df_plot["pme_TPM_y"] > 0]
df_plot["TPM_similarity"] = -abs(np.log(df_plot["pme_TPM_x"] / df_plot["pme_TPM_y"]))
df_plot["variance"] = df_plot.groupby("GeneID1")["pme_TPM_x"].transform("var")









df_plot
# Plot the distribution of pme_TPM_x vs pme_TPM_y, scatter plot
sns.scatterplot(data=df_plot, x="pme_TPM_x", y="pme_TPM_y", hue="Sequence Identity", s = 3)
# put same y and x limits
max_lim = max(df_plot["pme_TPM_x"].max(), df_plot["pme_TPM_y"].max())
plt.xlim(0, max_lim)
plt.ylim(0, max_lim)
plt.show()


# scatter plot of the sequence identity vs the expression similarity
sns.scatterplot(data=df_plot, x="Sequence Identity", y="TPM_similarity", hue= "pme_TPM_x", s=4, palette="viridis")
plt.show()


gene1 = "ENSG00000162998.4"

# --------------------------------------------------
# PROMOTERS WITH SEQUENCE IDENTITY == 1
# --------------------------------------------------
# select the promoters with sequence identity == 1 FROM DF PLOT 


genes = df_plot["GeneID1"].unique()
genes = np.random.choice(genes, 200)
# select the genes
test = df_plot[df_plot["GeneID1"].isin(genes)]
test = test[test["complete_name"].str.contains("liver")]
sns.scatterplot(data=test, x="Sequence Identity", y="TPM_similarity", hue= "pme_TPM_x", s=2, palette="viridis")
plt.show()


# make heatmap with all the variables in the dataframe and add star if it is significant
columns = ["Sequence Identity", "pme_TPM_x", "pme_TPM_y", "TPM_similarity"]
sns.heatmap(test[columns].corr(), annot=True)
# i cannot see the x labels
plt.xticks(rotation=45)
plt.show()

# compute the correlation between the sequence identity and the expression similarity
test[["Sequence Identity", "TPM_similarity"]].corr()


# make a boxplot of the expression similarity in bins of sequence identity
# bin the sequence identity
df_plot["Sequence Identity Bin"] = pd.cut(df_plot["Sequence Identity"], bins=10)
sns.boxplot(data=df_plot, x="Sequence Identity Bin", y="TPM_similarity")
plt.show()




# clean plot
plt.clf()
sns.scatterplot(data=df_plot, x="Sequence Identity", y="TPM_similarity", hue= "complete_name")
# y label 
plt.ylabel("Downstream gene expression similarity \n -abs(log(TPM similarity ratio))")
# put legend outside
plt.legend(bbox_to_anchor=(1.05, 1), loc='bottom right')
# change legend title
# Add title with the gene name 
plt.title("One promoter vs 70 other promoters")
plt.legend(title="sample")
plt.show()
# store plot
plt.savefig("./data/plots/sequence_identity_vs_expression_similarity_full.png")