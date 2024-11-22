

 python residPCA.py Initialize \
     --count_matrix_path /Users/shayecarver/residPCA/examples/example_data.h5ad \
     --vars_to_regress Batch,celltype,total_counts,pct_counts_mt,Age,Sex \
     --object_columns Batch,celltype,Sex \
     --variable_genes_flavor seurat \
     --n_PCs 150 \
     --random_seed 7 \
     --vargenes_IterPCA 3000 \
     --vargenes_Stand_resid 3000 \
     --BIC \
     --save_image_outputs \
     --basename test_run \
     --global_ct_cutoff 0.2
python residPCA.py Normalize --basename test_run
python residPCA.py Standardize --basename test_run
python residPCA.py StandardPCA_fit --basename test_run
