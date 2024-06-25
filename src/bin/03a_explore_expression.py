import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import seaborn as sns
import glob

# List with protein coding genes ids
pc_genes = pd.read_csv("./data/annotation/protein_coding_genes.gtf", sep="\t")
pc_genes_ids = pc_genes["gene_id"]

# get all files in the folder with absolute path
# use glob to get all files in the folder
expression_files = glob.glob("./data/expression/*/*.tsv")

# READ IN THE EXPERSSIONS, APPEND AND PLOT AND COLOR BY FILE NAME
exps = pd.DataFrame()
for file in expression_files:
    expression = pd.read_csv(file, sep="\t")
    name = file.split(".tsv")[0].split("/")[-1]
    print(name)
    expression["name"] = name
    expression["tissue"] = file.split("/")[-2]
    expression["complete_name"] = file.split("/")[-2] + "_" + name
    exps = pd.concat([exps, expression], ignore_index=True)

# keep only the protein coding genes
sample_pc = pc_genes_ids
exps = exps[exps["gene_id"].isin(sample_pc)]
# add pseudocount
exps["TPM"] = exps["TPM"] + 1
# plot histogram and log scale and color by name
# Get the seaborn color palette
blue_palette = sns.color_palette("Blues", n_colors=10)
green_palette = sns.color_palette("Reds", n_colors=10)
# Create a custom palette with the first two shades of blue and the second two shades of green
custom_palette = [blue_palette[6], blue_palette[9], green_palette[6], green_palette[9]]
# Set the palette
sns.set_palette(custom_palette)
sns.kdeplot(data=exps, x="TPM", hue="complete_name", log_scale=True)
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





