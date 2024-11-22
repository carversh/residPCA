#!/usr/bin/env python3

scExp = residPCA(
    count_matrix_path="./examples/example_data.h5ad",
    vars_to_regress=['Batch', 'celltype', 'total_counts', 'pct_counts_mt', 'Age', 'Sex'],
    object_columns=['Batch', 'Sex', 'celltype'],
    variable_genes_flavor="seurat",
    n_PCs=150,
    vargenes_IterPCA=3000,
    vargenes_Stand_resid=3000,
    BIC=True,
    save_image_outputs=True,
)

scExp.Normalize()

scExp.Standardize()

scExp.StandardPCA_fit()

scExp.ResidPCA_fit()

scExp.Iter_PCA_fit()