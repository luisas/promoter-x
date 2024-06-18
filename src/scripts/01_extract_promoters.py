# Script in python 
import pandas as pd
from pyfaidx import Fasta

gencode_gtf_path = "./data/annotation/gencode.v29.annotation.gtf"
mapper_entrez = "./data/mapper/gencode.v46.metadata.EntrezGene.gz"
fasta = Fasta('./data/genomes/GRCh38.primary_assembly.genome.fa')

# I took gencode annotation v29, the one most encode entries use. 
# Selected the entries that are protein coding and have a mapping with entrez id
# Took promoter regions as defined by 1kb upstream of the gene start and 0 downstream of the gene start


# --------------------------------------------------------------------
#                Read in files
# --------------------------------------------------------------------

# read in mapper as dictionary
mapper = pd.read_csv(mapper_entrez, sep="\t", header=None, names=["Ensembl", "EntrezID"])
mapper = dict(zip(mapper["Ensembl"], mapper["EntrezID"]))

# read in annotation
gtf_original = pd.read_csv(gencode_gtf_path, sep="\t", comment="#", header=None)
gtf_original["gene_id"] = gtf_original[8].str.extract(r'gene_id "(.*?)";')
gtf_original["transcript_id"] = gtf_original[8].str.extract(r'transcript_id "(.*?)";')
all_genes = gtf_original[gtf_original[2] == "gene"]
gtf_original = gtf_original[gtf_original[8].str.contains("protein_coding")]

all_protein_coding_genes = gtf_original[gtf_original[2] == "gene"]
# save all protein coding genes ids
all_protein_coding_genes.to_csv("./data/annotation/protein_coding_genes.gtf", sep="\t", index=False)

# --------------------------------------------------------------------
# Extract mapping from transcript id to entrez id
gtf = gtf_original[gtf_original[2] == "transcript"]
gtf["entrez_id"] = gtf["transcript_id"].map(mapper)
# remove rows with transcript id that are nas
gtf = gtf.dropna(subset=["entrez_id"])
gtf = gtf.drop_duplicates(subset="gene_id")
# keep only gene_ids and entrez_ids
mapped_genes = gtf[["gene_id", "entrez_id"]]


# --------------------------------------------------------------------
# extract mapped genes from gtf_original
gtf = gtf_original[gtf_original[8].str.contains("gene_id")]
# keep only genes
gtf = gtf[gtf[2] == "gene"]
# only keep genes that are in the mapped_genes
gtf = gtf[gtf["gene_id"].isin(mapped_genes["gene_id"])]
# keep columns chr, start, end, strand, gene_id, entrez_id
gtf = gtf[[0, 3, 4, 6, "gene_id"]]
# add entrez_id
gtf["entrez_id"] = gtf["gene_id"].map(dict(zip(mapped_genes["gene_id"], mapped_genes["entrez_id"])))
# rename columns
gtf.columns = ["Chromosome", "Gene Start", "Gene End", "Strand", "GeneID", "EntrezID"]
# make enrtez_id an integer
gtf["EntrezID"] = gtf["EntrezID"].astype(str)
# --------------------------------------------------------------------
# Extract promoter regions

upstream_threshold = 1000
downstream_threshold = 0 

# extract promoter regions and add to gtf
gtf["Promoter Start"] = gtf["Gene Start"] - upstream_threshold
gtf["Promoter End"] = gtf["Gene Start"] + downstream_threshold
gtf = gtf[gtf["Promoter Start"] > 0]
gtf = gtf[gtf["Promoter End"] > 0]


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
    promoters.drop(indices_to_drop, inplace=True)

# Call the function
remove_overlapping_promoters(gtf, all_genes)

# save these extracted promoters 
gtf.to_csv("./data/promoters/promoters.gtf", sep="\t", index=False)

gtf

# Create a fasta file with promoter regions
# extract sequence of promoter regions
with open('./data/promoters/promoters.fa', 'w') as f:
    # Iterate over each row in the DataFrame
    for i, row in gtf.iterrows():
        # Extract the sequence from the fasta file
        sequence = fasta[row["Chromosome"]][row["Promoter Start"]-1:row["Promoter End"]].seq
        # Write the sequence to the new fasta file
        f.write(f'>{row["GeneID"]}\n{sequence}\n')


# Use cd-hit to cluster regions and extract N representative promoters and extract M random promoters from each cluster
# This way we have a set of non-redundant promoters but also similar promoters (likely paralogs)

# Check the levels of gene expression of the downstream gene for each promoter (maybe it should be a pre-filtering step)



# Map pairs of promoters
# I want to extract N pairs of promoters that: 
# Span various levels of downstream gene expression 
# Span various levels of sequence similarity



# i first need to find non reduntant promoters
