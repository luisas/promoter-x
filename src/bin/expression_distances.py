#!/usr/bin/env python3

import numpy as np
from scipy.stats import pearsonr, spearmanr


def log_geometric_mean(sample_tpms): 
    log_tpms = np.log(sample_tpms)
    log_mean = np.mean(log_tpms)
    return log_mean

def pearson_correlation(x,y, log_sample_mean = None): 
    if log_sample_mean is not None: 
        x = np.log(x) - log_sample_mean
        y = np.log(y) - log_sample_mean
    return pearsonr(x,y)[0]

def spearman_correlation(x,y, log_sample_mean = None): 
    if log_sample_mean is not None: 
        x = np.log(x) - log_sample_mean
        y = np.log(y) - log_sample_mean
    return spearmanr(x,y)[0]

def corrected_logratio_mean(g1,g2,log_sample_mean):
    """
    Compute the corrected logratio mean of g1 and g2
    g1: gene expression of gene 1 across samples
    g2: gene expression of gene 2 across samples
    log_sample_mean: a metric of summary of gene expression for each sample
    """
    corrected_g1 = np.log(g1) - log_sample_mean
    mean_g1 = np.mean(corrected_g1)
    corrected_g2 = np.log(g2) - log_sample_mean
    mean_g2 = np.mean(corrected_g2)
    return abs(mean_g1 - mean_g2)



