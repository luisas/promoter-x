# Script in python 
import pandas as pd
from pyfaidx import Fasta
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

promoters = pd.read_csv("./data/promoters/promoters.gtf", sep="\t")
promoters_fasta = Fasta("./data/promoters/promoters.fa")

# extract 100 random promoters
promoters_one = promoters.sample(1)

# compute pairwise sequence identity between all promoters
promoter_pairs = []
for i, promoter1 in promoters_one.iterrows():
    for j, promoter2 in promoters.iterrows():
        if i < j:
            seq1 = str(promoters_fasta[str(promoter1["GeneID"])])
            seq2 = str(promoters_fasta[str(promoter2["GeneID"])])
            promoter_pairs.append([promoter1["GeneID"], promoter2["GeneID"], seq1, seq2])

promoter_pairs = pd.DataFrame(promoter_pairs, columns=["GeneID1", "GeneID2", "Sequence1", "Sequence2"])
promoter_pairs


# compute sequence identity with alignment
def compute_sequence_identity(seq1, seq2):
    alignments = pairwise2.align.globalxx(seq1, seq2)
    best_alignment = alignments[0]
    identity = best_alignment[2] / len(seq1)
    return identity


promoter_pairs["Sequence Identity"] = promoter_pairs.apply(lambda x: compute_sequence_identity(x["Sequence1"], x["Sequence2"]), axis=1)

promoter_pairs
# save it to a file
promoter_pairs.to_csv("./data/promoters/promoter_pairs.csv", index=False)
