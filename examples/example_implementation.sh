#!/bin/bash

source activate ResidPCA_Toolkit

# run from inside example folder of the github repo
residPCA Initialize \
     --count_matrix_path example_data.h5ad \
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

echo initialized

residPCA Normalize --basename test_run --path_to_directory ./
echo normalized 

residPCA Standardize --basename test_run --path_to_directory ./
echo standardized 

residPCA StandardPCA_fit --basename test_run --path_to_directory ./
echo performed standard PCA

residPCA residPCA_fit --basename test_run --path_to_directory ./
echo performed residPCA

residPCA Iter_PCA_fit --basename test_run --path_to_directory ./
echo performed Iterative PCA
residPCA ID_Global_CellType_States --basename test_run --path_to_directory ./

residPCA heritability_bed_output --basename test_run --path_to_directory ./ --ref_annotations_file ~/residPCA/gencode.v39.basic.annotation.names.bed --method Resid
echo outputtted heritability bed files from residPCA

residPCA heritability_bed_output --basename test_run --path_to_directory ./ --ref_annotations_file ~/residPCA/gencode.v39.basic.annotation.names.bed --method Standard
echo outputtted heritability bed files from  Standard PCA

residPCA heritability_bed_output --basename test_run --path_to_directory ./ --ref_annotations_file ~/residPCA/gencode.v39.basic.annotation.names.bed --method Iter
