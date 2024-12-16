#!/bin/bash

source activate ResidPCA_Toolkit

ResidPCA Initialize \
     --count_matrix_path ./examples/example_data.h5ad \
     --vars_to_regress Batch,celltype,total_counts,pct_counts_mt,Age,Sex \
     --object_columns Batch,celltype,Sex \
     --variable_genes_flavor seurat \
     --n_PCs 150 \
     --vargenes_IterPCA 3000 \
     --vargenes_Stand_resid 3000 \
     --BIC \
     --save_image_outputs \
     --path_to_directory ./ \
    --basename test_run


ResidPCA Normalize --basename test_run --path_to_directory ./

ResidPCA Standardize --basename test_run --path_to_directory ./

ResidPCA StandardPCA_fit --basename test_run --path_to_directory ./

ResidPCA residPCA_fit --basename test_run --path_to_directory ./

ResidPCA Iter_PCA_fit --basename test_run --path_to_directory ./

ResidPCA ID_Global_CellType_States --basename test_run --path_to_directory ./

ResidPCA heritability_bed_output --basename test_run --path_to_directory ./ --ref_annotations_file ~/residPCA/gencode.v39.basic.annotation.names.bed --method Resid

ResidPCA heritability_bed_output --basename test_run --path_to_directory ./ --ref_annotations_file ~/residPCA/gencode.v39.basic.annotation.names.bed --method Standard

ResidPCA heritability_bed_output --basename test_run --path_to_directory ./ --ref_annotations_file ~/residPCA/gencode.v39.basic.annotation.names.bed --method Iter