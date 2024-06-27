# Script in python 
import pandas as pd
from pyfaidx import Fasta
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import os

#               Read in promoter files 
# --------------------------------------------------
promoters_path = "./data/promoters/promoters.fa"
promoters_fasta = Fasta(promoters_path)
promoters = pd.read_csv("./data/promoters/promoters.gtf", sep="\t")

promoters

# --------------------------------------------------
# GET HIGHLY SIMILAR PROMOTERS
# --------------------------------------------------
prom_1 = "ENSG00000162998.4"
promoters_one = promoters[promoters["GeneID"] == prom_1]
path_promoters_sample_1 = "./data/promoters/sampled/promoters_sample_1.fa"
# store the sequences of the promoters in a fasta file
with open(path_promoters_sample_1, "w") as f:
    for i, promoter in promoters_one.iterrows():
        f.write(">"+str(promoter["GeneID"])+"\n")
        f.write(str(promoters_fasta[str(promoter["GeneID"])])+"\n")


# execute mmseqs command 
mmseqs_search_output = "./data/promoters/promoters_clustering_FULL.m8"
exec = "mmseqs easy-search "+ promoters_path +" "+ promoters_path+" "+mmseqs_search_output+" tmp --search-type 3"
os.system(exec)

# parse the output and take the first 100 promoters
similar_promoters = pd.read_csv(mmseqs_search_output, sep="\t", header=None)

# coun same entry in first column
prom_counts = similar_promoters[0].value_counts()
# store in a file
prom_counts.to_csv("./data/promoters/promoter_counts.csv", header=False)
seq_id = pd.read_csv(mmseqs_search_output, sep="\t", header=None)
seq_id = seq_id[[0, 1, 2]]
seq_id.columns = ["GeneID1", "GeneID2", "Sequence Identity"]
seq_id.to_csv("./data/promoters/promoter_pairs_test_FULL.csv", index=False)
seq_id

# --------------------------------------------------
# EXTRACT (POSSIBLY) DISSIMILAR PROMOTERS - SAMPLE RANDOMLY N
# --------------------------------------------------
# promoters
# random_promoters = promoters.sample(100)
# # append 
# promoters_ids = list(random_promoters["GeneID"]) + list(promoters_ids)
# promoters_two = promoters[promoters["GeneID"].isin(promoters_ids)]


# promoters_one
# promoters_two

# # search in with mmseqs one vs two
# path_promoters_sample_2 = "./data/promoters/sampled/promoters_sample_2.fa"
# # store the sequences of the promoters in a fasta file
# with open(path_promoters_sample_2, "w") as f:
#     for i, promoter in promoters_two.iterrows():
#         f.write(">"+str(promoter["GeneID"])+"\n")
#         f.write(str(promoters_fasta[str(promoter["GeneID"])])+"\n")

# # execute mmseqs command
# mmseqs_search_output = "./data/promoters/promoters_clustering_2.m8"
# exec = "mmseqs easy-search "+ path_promoters_sample_1 +" "+ path_promoters_sample_2 +" "+mmseqs_search_output+" tmp --search-type 3"
# os.system(exec)

# store in dataframe as GeneID1, GeneID2, Sequence identity (columns 0 1 2)

# # --------------------------------------------------                                                                                                                                                               
# # compute pairwise sequence identity between all promoters
# promoter_pairs = []
# for i, promoter1 in promoters_one.iterrows():
#     for j, promoter2 in promoters_two.iterrows():
#         if i < j:
#             seq1 = str(promoters_fasta[str(promoter1["GeneID"])])
#             seq2 = str(promoters_fasta[str(promoter2["GeneID"])])
#             promoter_pairs.append([promoter1["GeneID"], promoter2["GeneID"], seq1, seq2])

# promoter_pairs = pd.DataFrame(promoter_pairs, columns=["GeneID1", "GeneID2", "Sequence1", "Sequence2"])
# promoter_pairs



# # compute sequence identity with alignment
# def compute_sequence_identity(seq1, seq2):
#     alignments = pairwise2.align.globalxx(seq1, seq2)
#     best_alignment = alignments[0]
#     identity = best_alignment[2] / len(seq1)
#     return identity


# promoter_pairs["Sequence Identity"] = promoter_pairs.apply(lambda x: compute_sequence_identity(x["Sequence1"], x["Sequence2"]), axis=1)

# gene_test = "ENSG00000183307.3"
# promoter_pairs[(promoter_pairs["GeneID1"] == gene_test) | (promoter_pairs["GeneID2"] == gene_test)]
# promoter_pairs
# # save it to a file
# promoter_pairs.to_csv("./data/promoters/promoter_pairs_one.csv", index=False)
