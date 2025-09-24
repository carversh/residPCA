#!/usr/bin/env python3
import pandas as pd
import numpy as np
import scanpy
import scipy
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
import os
from datetime import datetime
from upsetplot import UpSet, plot, from_memberships
import pickle
from sklearn.decomposition import SparsePCA
import argparse
import re
# %%
import math
import statsmodels.api as sm


# %%
import sys
sys.path.append('/data/gusev/USERS/scarver/my_method')
from condPCA import *

def main(args):
    seed = args.seed
    method = args.method
    state_type = args.state_type
    dim = args.dim
    total_cells = args.total_cells
    perc_genes = args.perc_genes
    perc_cells = args.perc_cells
    flag_sim_100_cts = args.flag_sim_100_cts

    # Now you can use these parsed arguments as needed
    print("seed:", seed)
    print("method:", method)
    print("state_type:", state_type)
    print("dim:", dim)
    print("total_cells:", total_cells)
    print("perc_genes:", perc_genes)
    print("perc_cells:", perc_cells)
    print("flag_sim_100_cts:", flag_sim_100_cts)

# Creating an ArgumentParser instance
parser = argparse.ArgumentParser(description="Description of your program")

# Adding arguments
parser.add_argument("seed", type=int, help="Seed")
parser.add_argument("method", type=str, choices=["my_method","my_method_sparse"],
                    help="Method (scaled_NMF, PCA, Cond_PCA)")
parser.add_argument("state_type", type=str, choices=["within_one_ct", "across_cts"],
                    help="state_type (within_one_ct, across_cts)")
parser.add_argument("dim", type=int, help="Dimension (rank of matrix in DR method)")
parser.add_argument("total_cells",type=int, help="Number of cells (e.g., 'all' or a specific number)")
parser.add_argument("perc_genes",type=float, help="Percent genes in state")
parser.add_argument("perc_cells",type=float, help="Percent cells in state")
#parser.add_argument('flag_sim_100_cts', type=str, default='cts_7', help='flag 100 cts ')
parser.add_argument('--flag_sim_100_cts', type=str, default='cts_7', help='Optional flag_sim_100_cts')


# Parsing the arguments
args = parser.parse_args()

# Calling the main function with parsed arguments
main(args)



# %%
max_continuum = pd.read_csv('MAX_CONTINUUM_{}_{}_{}_{}_{}_gene_{}_cell_{}_flag_{}.txt'.format(
    args.method,
    args.total_cells,
    args.seed,
    args.state_type,
    args.dim,
    re.sub("[.]", "p", str(args.perc_genes)),
    re.sub("[.]", "p", str(args.perc_cells)),
    args.flag_sim_100_cts
), sep='\t', header=0, index_col=0)
# %%
ct_in_state = pd.read_csv('CT_IN_STATE_{}_{}_{}_{}_{}_gene_{}_cell_{}_flag_{}.txt'.format(
    args.method,
    args.total_cells,
    args.seed,
    args.state_type,
    args.dim,
    re.sub("[.]", "p", str(args.perc_genes)),
    re.sub("[.]", "p", str(args.perc_cells)),
    args.flag_sim_100_cts
), sep='\t', header=0, index_col=0).iloc[0,0]
# %%
metadata = pd.read_csv('METADATA_{}_{}_{}_{}_{}_gene_{}_cell_{}_flag_{}.txt'.format(
    args.method,
    args.total_cells,
    args.seed,
    args.state_type,
    args.dim,
    re.sub("[.]", "p", str(args.perc_genes)),
    re.sub("[.]", "p", str(args.perc_cells)),
    args.flag_sim_100_cts
), sep='\t', header=0, index_col=0)
# %%
my_method_light = condPCA(count_matrix_path='{}_{}_{}_{}_{}_gene_{}_cell_{}_flag_{}.txt'.format(
    args.method,
    args.total_cells,
    args.seed,
    args.state_type,
    args.dim,
    re.sub('[.]', 'p', str(args.perc_genes)),
    re.sub('[.]', 'p', str(args.perc_cells)),
    args.flag_sim_100_cts
),
metadata_path='METADATA_{}_{}_{}_{}_{}_gene_{}_cell_{}_flag_{}.txt'.format(
    args.method,
    args.total_cells,
    args.seed,
    args.state_type,
    args.dim,
    re.sub('[.]', 'p', str(args.perc_genes)),
    re.sub('[.]', 'p', str(args.perc_cells)),
    args.flag_sim_100_cts
),
object_columns=["celltype"],
vars_to_regress=["celltype"],
save_image_outputs=False,
BIC=False,
vargenes_Stand_Cond=2527,
n_PCs=args.dim,
random_seed=1,
vargenes_IterPCA=2527,
variable_genes_flavor="seurat",
global_ct_cutoff=0.2,
logged=True)

# %%
# log normalize count data
my_method_light.Normalize()
# # standardize count data and metadata
my_method_light.Standardize()
# # fit standard PCA
my_method_light.StandardPCA_fit()
my_method_light.CondPCA_fit()
#pd.DataFrame(my_method_light.standardized_residual).to_csv('my_method_cond_resid.csv', index=True)
#my_method_light.CondPCA_cell_embeddings.to_csv('my_method_cond_emb.csv', index=True)
# %%
def adj_rsq_compute(X_input):
    scaled_emb = (X_input - X_input.mean()) / X_input.std()
    scaled_time = (max_continuum["State_1"] - max_continuum["State_1"].mean()) / max_continuum["State_1"].std()
    # indices dont match so rename according to cell names (simulated)
    scaled_time.index = scaled_emb.index

    # Add constant column to the predictor variables
    X = sm.add_constant(scaled_emb)

    # Fit linear regression model
    model = sm.OLS(scaled_time, X).fit()

    # Extract the adjusted R-squared value from the summary
    adj_r_squared = model.rsquared_adj
    return adj_r_squared

def max_correlation_compute(X_input):
    scaled_emb = (X_input - X_input.mean()) / X_input.std()
    scaled_time = (max_continuum["State_1"] - max_continuum["State_1"].mean()) / max_continuum["State_1"].std()
    # indices dont match so rename according to cell names (simulated)
    scaled_time.index = scaled_emb.index
    # Initialize an empty list to store adjusted R-squared values
    corr = []

    # Fit separate linear regression models for each predictor variable
    for column in scaled_emb.columns:
        # Add constant column to predictor variable
        X = sm.add_constant(scaled_emb[column])
        # Fit linear regression model
        model = sm.OLS(scaled_time, X).fit()
        # Append adjusted R-squared value to list
        corr.append(model.rsquared_adj)
    corr_series = pd.Series(corr, index=scaled_emb.columns)
    return np.max(corr_series)

def max_correlation_compute_SUB(X_input):
    #print(X_input.shape)
    sub_embeddings = X_input.loc[(metadata["celltype"] == ct_in_state).values,:]
    #sub_embeddings = X_input.loc[metadata["celltype"] == ct_in_state,:]
    sub_state = max_continuum.loc[(metadata["celltype"] == ct_in_state).values,:]
    # indices dont match so rename according to cell names (simulated)
    sub_state.index = sub_embeddings.index
    # Initialize an empty list to store adjusted R-squared values
    corr = []

   # Fit separate linear regression models for each predictor variable
    for column in sub_embeddings.columns:
        # Add constant column to predictor variable
        correlation = np.corrcoef(sub_embeddings[column], sub_state["State_1"])[0, 1]
        # Append adjusted R-squared value to list
        corr.append(correlation**2)
    corr_series = pd.Series(corr, index=sub_embeddings.columns)
    return np.max(corr_series)
        # Add constant column to predictor variable
   #     X = sm.add_constant(sub_embeddings[column])
        # Fit linear regression model
    #    model = sm.OLS(sub_state, sub_embeddings).fit()
        # Append adjusted R-squared value to list
     #   corr.append(model.rsquared_adj)
   # corr_series = pd.Series(corr, index=scaled_emb.columns)
   # return np.max(corr_series)

# %%
print(adj_rsq_compute(my_method_light.CondPCA_cell_embeddings))
print(adj_rsq_compute(my_method_light.StandardPCA_cell_embeddings))
cond_adj = adj_rsq_compute(my_method_light.CondPCA_cell_embeddings)
stand_adj = adj_rsq_compute(my_method_light.StandardPCA_cell_embeddings)

print(max_correlation_compute(my_method_light.CondPCA_cell_embeddings))
print(max_correlation_compute(my_method_light.StandardPCA_cell_embeddings))
cond_max_corr = max_correlation_compute(my_method_light.CondPCA_cell_embeddings)
stand_max_corr = max_correlation_compute(my_method_light.StandardPCA_cell_embeddings)

#print(my_method_light.CondPCA_cell_embeddings.shape)
#print(my_method_light.CondPCA_gene_loadings.shape)

if args.state_type == "within_one_ct":
    SUB_cond_max_corr = max_correlation_compute_SUB(my_method_light.CondPCA_cell_embeddings)
    SUB_stand_max_corr = max_correlation_compute_SUB(my_method_light.StandardPCA_cell_embeddings)
else:
    SUB_cond_max_corr = cond_max_corr
    SUB_stand_max_corr = stand_max_corr

# %%

data = {
    "seed": [args.seed],
    "method": ["my_conPCA"],
    "state_type": [args.state_type],
    "dim": [args.dim],
    "total_cells": [args.total_cells],
    "perc_genes": [args.perc_genes],
    "perc_cells": [args.perc_cells],
    "flag": [args.flag_sim_100_cts],
    "ct_in_state": [ct_in_state],
    "adj.rsq": [cond_adj],
    "max.rsq": [cond_max_corr],
    "max.rsq_sub": [SUB_cond_max_corr],
    "total_cells_can_occupy_state": [math.floor(my_method_light.StandardPCA_cell_embeddings.loc[(metadata["celltype"] == ct_in_state).values,:].shape[0] )] if args.state_type == "within_one_ct" else [my_method_light.StandardPCA_cell_embeddings.shape[0]],
    "total_cells_occupy_state": [math.floor(my_method_light.StandardPCA_cell_embeddings.loc[(metadata["celltype"] == ct_in_state).values,:].shape[0] * args.perc_cells)] if args.state_type == "within_one_ct" else [my_method_light.StandardPCA_cell_embeddings.shape[0]* args.perc_cells],
    "total_genes": [my_method_light.StandardPCA_gene_loadings.shape[0]],
    "total_genes_in_state": [math.floor(my_method_light.StandardPCA_gene_loadings.shape[0] * args.perc_genes)]
}



# %%
# Create DataFrame from dictionary
df = pd.DataFrame(data)

# %%
# Write DataFrame to file, append mode
with open("output_my_method.txt", "a") as f:
    df.to_csv(f, sep=" ", header=False, index=False)

# %%

data = {
    "seed": [args.seed],
    "method": [args.method],
    "state_type": [args.state_type],
    "dim": [args.dim],
    "total_cells": [args.total_cells],
    "perc_genes": [args.perc_genes],
    "perc_cells": [args.perc_cells],
    "flag": [args.flag_sim_100_cts],
    "ct_in_state": [ct_in_state],
    "adj.rsq": [stand_adj],
    "max.rsq": [stand_max_corr],
    "max.rsq_sub": [SUB_stand_max_corr],
    "total_cells_can_occupy_state": [math.floor(my_method_light.StandardPCA_cell_embeddings.loc[(metadata["celltype"] == ct_in_state).values,:].shape[0] )] if args.state_type == "within_one_ct" else [my_method_light.StandardPCA_cell_embeddings.shape[0]],
    "total_cells_occupy_state": [math.floor(my_method_light.StandardPCA_cell_embeddings.loc[(metadata["celltype"] == ct_in_state).values,:].shape[0] * args.perc_cells)] if args.state_type == "within_one_ct" else [my_method_light.StandardPCA_cell_embeddings.shape[0]* args.perc_cells],
    "total_genes": [my_method_light.StandardPCA_gene_loadings.shape[0]],
    "total_genes_in_state": [math.floor(my_method_light.StandardPCA_gene_loadings.shape[0] * args.perc_genes)]
}

# %%
# Create DataFrame from dictionary
df = pd.DataFrame(data)

# Write DataFrame to file, append mode
with open("output_my_method.txt", "a") as f:
    df.to_csv(f, sep=" ", header=False, index=False)


# deleting tmp files
os.remove('MAX_CONTINUUM_{}_{}_{}_{}_{}_gene_{}_cell_{}_flag_{}.txt'.format(
    args.method,
    args.total_cells,
    args.seed,
    args.state_type,
    args.dim,
    re.sub("[.]", "p", str(args.perc_genes)),
    re.sub("[.]", "p", str(args.perc_cells)),
    args.flag_sim_100_cts
))
os.remove('CT_IN_STATE_{}_{}_{}_{}_{}_gene_{}_cell_{}_flag_{}.txt'.format(
    args.method,
    args.total_cells,
    args.seed,
    args.state_type,
    args.dim,
    re.sub("[.]", "p", str(args.perc_genes)),
    re.sub("[.]", "p", str(args.perc_cells)),
    args.flag_sim_100_cts
))
os.remove('METADATA_{}_{}_{}_{}_{}_gene_{}_cell_{}_flag_{}.txt'.format(
    args.method,
    args.total_cells,
    args.seed,
    args.state_type,
    args.dim,
    re.sub("[.]", "p", str(args.perc_genes)),
    re.sub("[.]", "p", str(args.perc_cells)),
    args.flag_sim_100_cts
))
os.remove('{}_{}_{}_{}_{}_gene_{}_cell_{}_flag_{}.txt'.format(
    args.method,
    args.total_cells,
    args.seed,
    args.state_type,
    args.dim,
    re.sub('[.]', 'p', str(args.perc_genes)),
    re.sub('[.]', 'p', str(args.perc_cells)),
    args.flag_sim_100_cts
))
