# Script in python 
import pandas as pd
from pyfaidx import Fasta
import matplotlib.pyplot as plt

gencode_gtf_path = "./data/annotation/gencode.v29.annotation.gtf"
mapper_entrez = "./data/mapper/gencode.v46.metadata.EntrezGene.gz"
fasta = Fasta('./data/genomes/GRCh38.primary_assembly.genome.fa')

# I took gencode annotation v29, the one most encode entries use. 
# Selected the entries that are protein coding and have a mapping with entrez id
# Took promoter regions as defined by 1kb upstream of the gene start and 0 downstream of the gene start


# --------------------------------------------------------------------
#                Read in files
# --------------------------------------------------------------------

# read in annotation
gtf_original = pd.read_csv(gencode_gtf_path, sep="\t", comment="#", header=None)
gtf_original["gene_id"] = gtf_original[8].str.extract(r'gene_id "(.*?)";')
gtf_original = gtf_original[gtf_original[8].str.contains("protein_coding")]
all_protein_coding_genes = gtf_original[gtf_original[2] == "gene"]
# save all protein coding genes ids
all_protein_coding_genes.to_csv("./data/annotation/protein_coding_genes.gtf", sep="\t", index=False)


gtf_original = gtf_original[gtf_original[2] == "gene"]
gtf = gtf_original[[0, 3, 4, 6, "gene_id"]]
gtf.columns = ["Chromosome", "Gene Start", "Gene End", "Strand", "GeneID"]

gtf

def remove_overlapping_promoters(promoters, all_genes):
    indices_to_drop = []
    # Iterate over each row in the DataFrame
    for i, promoter in promoters.iterrows():
        # Check if the promoter overlaps with any gene
        same_chromosome = all_genes[0] == promoter['Chromosome']
        overlap_start = all_genes[3] <= promoter['Promoter End']
        overlap_end = all_genes[4] >= promoter['Promoter Start']
        overlaps = same_chromosome & overlap_start & overlap_end

        # If the promoter overlaps with any gene other than its own, add its index to the list
        if any(overlaps) and len(all_genes[overlaps]['gene_id'].unique()) > 1:
            indices_to_drop.append(i)

    # Drop the overlapping promoters from the DataFrame
    promoters = promoters.drop(indices_to_drop)
    return(promoters)

def get_promoters(gtf, upstream_threshold, downstream_threshold):

    # extract promoter regions and add to gtf
    gtf["Promoter Start"] = gtf["Gene Start"] - upstream_threshold
    gtf["Promoter End"] = gtf["Gene Start"] + downstream_threshold
    gtf = gtf[gtf["Promoter Start"] > 0]
    gtf = gtf[gtf["Promoter End"] > 0]

    gtf = remove_overlapping_promoters(gtf, all_protein_coding_genes)
    return(gtf)

def count_promoters(gtf):
    return(gtf.shape[0])


# ------------------------------------------------
#               Extract promoters
# ------------------------------------------------
upstream_threshold = 100
downstream_threshold = 0 
gtf = get_promoters(gtf, upstream_threshold, downstream_threshold)




# ------------------------------------------------
# Save and extract
# ------------------------------------------------

# save these extracted promoters 
gtf.to_csv("./data/promoters/promoters_100.gtf", sep="\t", index=False)

# Create a fasta file with promoter regions
# extract sequence of promoter regions
with open('./data/promoters/promoters_100.fa', 'w') as f:
    # Iterate over each row in the DataFrame
    for i, row in gtf.iterrows():
        # Extract the sequence from the fasta file
        sequence = fasta[row["Chromosome"]][row["Promoter Start"]-1:row["Promoter End"]].seq
        # Write the sequence to the new fasta file
        f.write(f'>{row["GeneID"]}\n{sequence}\n')



# ---------------------------------------------------------------------------
# Benchmark size of promoter effect on number of promoters
# this is the consequence of the overlap filter 
# where we remove promoters that overlap with other protein coding genes
# ---------------------------------------------------------------------------
# check the number of promoters vs size of threshold
thresholds = [1000, 500, 100, 50, 5]
promoter_count_benchmark = []
for threshold in thresholds:
    print(threshold)
    promoter_count_benchmark.append(count_promoters(get_promoters(gtf, threshold, 0)))

# save the counts in a file
promoter_count_benchmark_df = pd.DataFrame({"Threshold": thresholds, "Promoter Count": promoter_count_benchmark})
promoter_count_benchmark_df.to_csv("./data/promoters/promoter_count_benchmark.csv", index=False)

promoter_count_benchmark_df

#  plot the number of promoters vs size of threshold
plt.plot(thresholds, promoter_count_benchmark)
# scale limit y from 0 to 16000
plt.ylim(0, 16000)
plt.xlabel("Threshold")
plt.ylabel("Number of promoters")
plt.show()