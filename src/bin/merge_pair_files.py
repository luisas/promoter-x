#!/usr/bin/env python3

import pandas as pd
import sys


file1 = sys.argv[1]
file2 = sys.argv[2]
outfile = sys.argv[3]

df1 = pd.read_csv(file1)
df2 = pd.read_csv(file2)

# get gene1 and gene2, sort lexicographically btw the two, and paste them together, then add a new column 

def create_id(row):
    gene1 = row['gene1']
    gene2 = row['gene2']
    if gene1 < gene2:
        return gene1 + '-' + gene2
    else:
        return gene2 + '-' + gene1

df1['pair_id'] = df1.apply(create_id, axis=1)
df1 = df1.drop_duplicates(subset='pair_id')
df2['pair_id'] = df2.apply(create_id, axis=1)
df2 = df2.drop_duplicates(subset='pair_id')
# remove the gene1 and gene2 columns with _y suffix
df2 = df2.drop(columns=['gene1', 'gene2'])

# merge the two dataframes on the pair_id column, outer join
merged_df = pd.merge(df1, df2, on='pair_id', how='outer')

# save the merged dataframe to a file
merged_df.to_csv(outfile, index=False)
