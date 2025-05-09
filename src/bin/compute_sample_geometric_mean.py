#!/usr/bin/env python3
import os
import pandas as pd
import glob
import sys
import numpy as np
np.random.seed(42)
import expression_distances as dist

expression_dir = sys.argv[1]
expression_files = glob.glob(expression_dir+"/**/*.tsv")
outfile = sys.argv[2]


expression_dir = "/home/luisasantus/Desktop/crg_cluster/projects/PROMOTER-X/data/expression"

log_geom_means = []
names = []
for file in expression_files:
    expression = pd.read_csv(file, sep="\t")
    name = file.split(".tsv")[0].split("/")[-1]
    complete_name = file.split("/")[-2] + "_" + name
    pme_TPM = expression["pme_TPM"]+1
    names.append(complete_name)
    log_geom_means.append(dist.log_geometric_mean(pme_TPM))

# save the log geometric means to a csv file with name
log_geom_means = pd.DataFrame(log_geom_means, columns=["log_geom_mean"])
log_geom_means["name"] = names
log_geom_means.to_csv(outfile, index = False)