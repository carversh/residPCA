#!/usr/bin/env python3
# ensure that you have activated your conda environment before running this script.
# example:
# source activate ResidPCA_Toolkit
# python example_implementation.py

from residPCA import residPCA

scExp = residPCA(
    count_matrix_path="example_data.h5ad",
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

scExp.residPCA_fit()

scExp.Iter_PCA_fit()

scExp.heritability_bed_output(ref_annotations_file = "~/residPCA/gencode.v39.basic.annotation.names.bed", num_genes = 200, method = "Resid")

scExp.heritability_bed_output(ref_annotations_file = "~/residPCA/gencode.v39.basic.annotation.names.bed", num_genes = 200, method = "Standard")

scExp.heritability_bed_output(ref_annotations_file = "~/residPCA/gencode.v39.basic.annotation.names.bed", num_genes = 200, method = "Iter")

# access cell embeddings
scExp.residPCA_cell_embeddings
scExp.StandardPCA_cell_embeddings
scExp.IterPCA_cell_embeddings
# access gene loadings
scExp.residPCA_gene_loadings
scExp.StandardPCA_gene_loadings
scExp.IterPCA_gene_loadings
# access BIC cutoffs
scExp.residPCA_BIC_cutoff
scExp.StandardPCA_BIC_cutoff
scExp.IterPCA_BIC_cutoff
