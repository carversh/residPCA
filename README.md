# ResidPCA Package v1.0.0

The Residual Principal Component Analysis (ResidPCA) Toolkit is a comprehensive software platform designed to implement a novel method, ResidPCA, and well as Standard PCA and Iterative PCA (PCA applied to each cell type) to identify sets of cellular states within single cell data. ResidPCA leverages log-normalized TP10k data with known cell type labels to estimate a set of cellular states that are independent of cell type. The method first removes cell type-driven noise from the expression matrix and then applies PCA on the residucalized matrix to identify a set of denoised cellular states within and across cell types.

As input, the pipline accepts a single-cell cell x genes count matrix (either a single-cell RNA-seq matrix or a single-cell ATAC-seq matrix with peaks mapped to genes) in the form of a [.h5ad](https://anndata.readthedocs.io/en/latest/tutorials/notebooks/getting-started.html) file or a .txt file. The accompanying metadata (including cell type labels) must be included in the .h5ad object or inputted as a separate .txt file. 

Note: standard pre-processing pipelines can map scATAC-seq data to peaks such as Signac package in Seurat.

The pre-print is available [here]().

![Inline Image](https://github.com/carversh/residPCA/blob/main/residPCA_visual.png)

# Installation

First, install the neccessary conda environment in which to run the ResidPCA Toolkit (the environment configuration file can be found [here](https://github.com/carversh/residPCA/blob/main/environment.yml) ):
```
conda env create -f environment.yml
```
Note: ensure that [Miniconda](https://docs.anaconda.com/miniconda/install/) or [Conda](https://anaconda.org/anaconda/conda) is installed within your computing system.

Second, activate the Conda environment:
```
source activate ResidPCA_Toolkit
```
 
This conda environment now contains all the necessary packages to run the Toolkit.

Third, install the ResidPCA Toolkit via [pip](https://pypi.org/):

```
pip install ResidPCA
```

Note: when running the ResidPCA Toolkit, the respective conda environment must always be activated using the command:
```
source activate ResidPCA_Toolkit
```
To deactivate this environment run:
```
conda deactivate ResidPCA_Toolkit
```
# Tutorials

Command Line Interface:
  - A bash script demonstrating an example run of the pipeline for **simulated data** is available [here]().
  - A bash script for **real data** is available [here](https://github.com/carversh/residPCA/blob/main/examples/example_implementation.sh).

Python Environment:
  - A Python script performing the same task for **simulated data** is provided [here]().
  - A Python script for **real data** is provided [here](https://github.com/carversh/residPCA/blob/main/examples/example_implementation.py).

The next section offers a detailed, step-by-step guide explaining each command in the process.

# Step by Step Guide 

The ResidPCA Toolkit can be run from the command line or the same steps can be run from within a Python environment. In this tutorial, we will show you how to run the method using either option.

### Step 1 - instantiate a class with your input data

Command line example  command:
```
ResidPCA Initialize \
     --count_matrix_path ./examples/example_data.h5ad \
     --vars_to_regress Batch,celltype,total_counts,pct_counts_mt,Age,Sex \
     --object_columns Batch,celltype,Sex \
     --variable_genes_flavor seurat \
     --n_PCs 150 \
     --vargenes_IterPCA 3000 \
     --vargenes_Stand_resid 3000 \
     --BIC \
     --save_image_outputs
```

Python environment example command:
```
scExp = residPCA(
    count_matrix_path="./examples/example_data.h5ad",
    vars_to_regress=['Batch', 'celltype', 'total_counts', 'pct_counts_mt', 'Age', 'Sex'],
    object_columns=['Batch', 'Sex', 'celltype'],
    variable_genes_flavor="seurat",
    n_PCs=150,
    vargenes_IterPCA=3000,
    vargenes_Stand_resid=3000,
    BIC=True,
    save_image_outputs=True
)
```

**Input Data**
  - ```count_matrix_path``` - points to the .h5ad object or tab delimited file that is cells by genes dimensional, where rownames correspond to the barcoded cells and columnnames correspond to the genes.
  - ```metadata_path``` - IF metadata is not included in the .h5ad object, this parameter points to the tab delimited file that is cells by covariates/features dimensional.There must be a single column with a column name "celltype" that contains the cell type labels corresponding to each barcoded cell in the count matrix.

**Parameters:**
- `--count_matrix_path` (str, required):
Path to the count matrix file.
`--metadata_path` (str, optional):
Path to the metadata file. Default is None.
- `--vars_to_regress` (str, optional):
Comma-separated list of variables to regress. Default is True.
- `--object_columns` (str, required):
Comma-separated list of object columns.
- `--variable_genes_flavor` (str, optional):
Specifies the flavor of variable genes to use. Options include:
    -`seurat`
    -`cell_ranger`
    -`seurat_v3`
    -`seurat_v3_paper`
Default is seurat_v3.
- `--n_PCs` (int, optional):
Number of principal components to compute. Default is 200.
- `--random_seed` (int, optional):
Random seed for reproducibility. Default is 9989999.
- `--vargenes_Stand_resid` (int or str, optional):
Number of variable genes to use for Standard PCA and ResidPCA. Accepts an integer or all. Default is all meaning all genes are included in analysis.
- `--vargenes_IterPCA` (int or str, optional):
Number of variable genes to use for Iterative PCA. Accepts an integer or all. Default is all meaning all genes are included in analysis.
- `--BIC` (bool, optional):
Use BIC for model selection. Enabled by default.
- `--no_BIC` (bool, optional):
Do not use BIC for model selection. Overrides --BIC.
- `--save_image_outputs` (bool, optional):
Save image outputs. Disabled by default.
- `--path_to_directory` (str, optional):
Path to the output directory. Default is "./".
- `--basename` (str, optional):
Basename for output files. Default is residPCA_run_<current_datetime>.
- `--global_ct_cutoff` (float, optional):
Global cutoff for cell types. Default is 0.2.
`--logged` (bool, optional):
Indicate if data is logged. Disabled by default.
- `--sparse_PCA` (bool, optional):
Use sparse PCA. Disabled by default.

### Step 2 - log-normalize the count data

Command line example  command:
```
ResidPCA Normalize 
```

Python environment example command:
```
scExp.Normalize()
```

**Parameters:**
- `--path_to_directory` (str, optional):
Path to the output directory. Default is "./".
- `--basename` (str, optional):
Basename for output files. Default is residPCA_run_<current_datetime>.

### Step 3 - standardize the count data

Command line example  command:
```
ResidPCA Standardize 
```
Python environment example command:
```
scExp.Standardize()
```
**Parameters:**
- `--path_to_directory` (str, optional):
Path to the output directory. Default is "./".
- `--basename` (str, optional):
Basename for output files. Default is residPCA_run_<current_datetime>.

### Step 4 - perform Standard PCA

Command line example  command:
```
ResidPCA StandardPCA_fit 
```

Python environment example command:
```
scExp.StandardPCA_fit()
```
Returns the Standard PCA output in the form of dataframes. 

**Parameters:**
- `--path_to_directory` (str, optional):
Path to the output directory. Default is "./".
- `--basename` (str, optional):
Basename for output files. Default is residPCA_run_<current_datetime>.

**Outputs in command line:**
  -`StandardPCA_gene_loadings.csv`
  -`StandardPCA_cell_embeddings.csv`
Note: embeddings/loadings output at BIC cuttoff or n_PCs specified

**Outputs in Python object:**
  - ```scExp.StandardPCA_cell_embeddings``` - cell embeddings outputted by Standard PCA
  - ```scExp.StandardPCA_gene_loadings``` - gene loadings or eigenvectors outputted by Standard PCA
  - ```scExp.StandardPCA_eigenvalues``` - eigenvalues outputted by Standard PCA
  - ```scExp.StandardPCA_BIC_cutoff``` - PC cutoff that specifies the maximum state that is significant. For significant states, subset the cell embeddings and gene loadings from PC1 to the PC specified in this variable

### Step 5 - perform Residual PCA (ResidPCA)

Command line example  command:
```
ResidPCA ResidPCA_fit 
```

Python environment example command:
```
scExp.ResidPCA_fit()
```
Returns the CondPCA output in the form of dataframes. 

**Parameters:**
- `--path_to_directory` (str, optional):
Path to the output directory. Default is "./".
- `--basename` (str, optional):
Basename for output files. Default is residPCA_run_<current_datetime>.

**Outputs in command line:**
  -`ResidPCA_gene_loadings.csv`
  -`ResidPCA_cell_embeddings.csv`
  Note: embeddings/loadings output at BIC cuttoff or n_PCs specified

**Outputs in Python object:**
  - ```scExp.CondPCA_cell_embeddings``` - cell embeddings outputted by Conditional PCA
  - ```scExp.CondPCA_gene_loadings``` - gene loadings or eigenvectors outputted by Conditional PCA
  - ```scExp.CondPCA_eigenvalues``` - eigenvalues outputted by Conditional PCA
  - ```scExp.CondPCA_BIC_cutoff``` - PC cutoff that specifies the maximum state that is significant. For significant states, subset the cell embeddings and gene loadings from PC1 to the PC specified in this variable

### Step 6 - perform Iterative PCA (IterPCA)

Command line example  command:
```
ResidPCA Iter_PCA_fit 
```

Python environment example command:
```
scExp.Iter_PCA_fit()
```
Returns the IterPCA output in the form of dictionaries, where each dictionary is equal to the length of the number of cell types. The keys of the dictionary correspond to the cell type in the "celltype" column of the metadata while the values of the dictionary represent the respective dataframe that corresponds to that cell type. 

**Outputs in command line** (the following files will be generated as outputs):

- `Iter_PCA_cell_embeddings_*.csv`: contains cell embeddings data for each iteration. 
- `Iter_PCA_gene_loadings_*.csv`: contains gene loadings data for each iteration.

The state embeddings and gene loadings are outputted at the BIC cutoff if the BIC cutoff is flagged. If not, the specified number of PCs parameterized by `n_PCs` will be used.
The * in the filenames represents the name of the cell type as specified in the metadata. For each cell type, a pair of files (gene loadings and cell embeddings) will be outputted.

**Outputs in Python object:**
  - ```scExp.IterPCA_cell_embeddings``` - dictionary containing the cell embeddings outputted by Iterative PCA per cell type
  - ```scExp.IterPCA_gene_loadings``` - dictionary containing the gene loadings outputted by Iterative PCA per cell type
  - ```scExp.IterPCA_eigenvalues``` - dictionary containing the gene eigenvalues outputted by Iterative PCA per cell type
  - ```scExp.IterPCA_BIC_cutoff``` - dictionary of PC cutoffs that specifies the maximum state that is significant per cell type. For significant states, subset the cell embeddings and gene loadings from PC1 to the PC specified in this variable

Warning:
  - If there are fewer than 200 cells in a cell type, the method will return an empty dataset for that given cell type.

### Step 7 - identify states that are cell type specific and states that span all cell types

Command line example  command:
```
ResidPCA ID_Global_CellType_States 
```

Python environment example command:
```
scExp.ID_Global_CellType_States()
```

Returns a dataframe for StandardPCA based states and CondPCA based states. Returns the maximum squared correlation between states identified in IterPCA and StandardPCA/CondPCA. Annotates each state and whether it is a global or cell type specific state as well as the cell types that are in that state.

Note: This method labels each state identified in StandardPCA and CondPCA as a global state, or a state spanning multiple cell types, or cell type specific states, or a state within a specific cell type. Additionally, this method labels which states belong to which cell types. To perform this method, ```scExp.Iter_PCA_fit()``` must be run beforehand and at least ```scExp.StandardPCA_fit()``` or ```scExp.CondPCA_fit()```  must be run. If only ```scExp.CondPCA_fit()``` is run, only the states belonoging to CondPCA will be evaluated, if only ```scExp.StandardPCA_fit()``` is run, only the states belonging to StandardPCA will be evaluated, if both are run, both will be evaluated.

**Outputs **
 -```scExp.StandardPCA_IterPCA_squared_correlations``` - dataframe containing each state from StandardPCA and a row and the maximum correlation with a certain run of IterPCA (columns are labeled by IterPCA run on that given cell type). Two additional columns are included called "CT_involved", which tells you which cell types are involved in the state, and "Global_vs_CT", which tell you whether the state is global or cell type specific.
 - ```scExp.CondPCA_IterPCA_squared_correlations``` - dataframe containing each state from CondPCA and a row and the maximum correlation with a certain run of IterPCA (columns are labeled by IterPCA run on that given cell type). Two additional columns are included called "CT_involved", which tells you which cell types are involved in the state, and "Global_vs_CT", which tell you whether the state is global or cell type specific.

 If ```save_image_outputs == True```:
   - example
   

### Step 9 - Output Loadings from Each Method as BED Files for LDSC-SEG (Heritability Analysis)

The loadings from each dimensionality reduction method can be extracted and saved as BED files. These BED files annotate the top (most positive) and bottom (most negative) gene loadings, which can then be used as input for **LDSC-SEG** (partitioned heritability analysis).

### Purpose
The exported BED files are used to annotate genomic regions corresponding to the top and bottom gene loadings, enabling functional genomic analyses and heritability partitioning.

Command line example  command:
```
ResidPCA heritability_bed_output --basename test_run --path_to_directory ./ --ref_annotations_file ~/residPCA/gencode.v39.basic.annotation.names.bed --method Resid # output bed files for loadings from ResidPCA

ResidPCA heritability_bed_output --basename test_run --path_to_directory ./ --ref_annotations_file ~/residPCA/gencode.v39.basic.annotation.names.bed --method Standard # output bed files for loadings from

ResidPCA heritability_bed_output --basename test_run --path_to_directory ./ --ref_annotations_file ~/residPCA/gencode.v39.basic.annotation.names.bed --method Iter # output bed files for loadings from IterPCA
```

Python environment example command:
```
scExp.heritability_bed_output("~/residPCA/gencode.v39.basic.annotation.names.bed", 200, "Resid") # output bed files for loadings from ResidPCA

scExp.heritability_bed_output("~/residPCA/gencode.v39.basic.annotation.names.bed", 200, "Standard") # output bed files for loadings from StandardPCA

scExp.heritability_bed_output("~/residPCA/gencode.v39.basic.annotation.names.bed", 200, "Iter") # output bed files for loadings from IterPCA
```

**Input Data**
  - ```ref_annotations_file``` - This is the reference file containing genomic annotations, such as gene names and chromosomal positions.
Example file: [gencode.v39.basic.annotation.names.bed](https://github.com/carversh/residPCA/blob/main/gencode.v39.basic.annotation.names.bed) (provided in the repository).
  It contains:
    Gene names.
    Chromosomal start and end positions.
    Additional metadata for annotations.

**Parameters:**
- `--ref_annotations_file` (str, required):
Path to the file containing genomic annotations. 
`--num_genes` (str, optional):
Specifies the number of genes to include in the annotation. Default: 200.
- `--method` (str, optional):
Specifies the dimensionality reduction method for which BED files should be generated. Options include:
    -`Resid` (for ResidPCA)
    -`Standard` (for Standard PCA)
    -`Iter` (for Iterative PCA)
- `--window` (str, optional):
Window size (in base pairs) for annotations upstream and downstream of the gene. Default: 100,000 (100 kb).


# Image Outputs

When initializing your experiment, if you set ```save_image_outputs = True```, a directory will be created that will contain all image based outputs from the method. Different commands will yield different image outputs depending on the task of the command. All relevant images and data will be saved to a directory called ```basename``` with the initial appended path ```path_to_directory```.

# Low Memory Setting


Here’s an improved version of your write-up, providing more clarity, structure, and a professional tone:

The low-memory mode is designed to reduce memory usage by limiting the number of intermediate outputs saved during computation. As a result, it can only run one method at a time. For instance, low-memory mode can execute either ResidPCA or Standard PCA in a single instance, but it cannot handle both simultaneously.

Currently, low-memory support for Iterative PCA is still under development.

In this mode, several intermediate matrices are not stored in memory to conserve resources. These include:

- Normalized count matrices
- Standardized count matrices
- Metadata
- The residualized matrix (specific to ResidPCA)

Despite these limitations, the key outputs—cell embeddings and gene loadings—are still saved and remain accessible for downstream analyses.

If you want to run ResidPCA, instantiate your class with the ```lowmem=True``` and either run ```ResidPCA_fit()``` or ```StandardPCA_fit()```.

Command line example  command:
```
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
     --basename test_run_LOWMEM \
     --path_to_directory ./ \
     --lowmem
     
ResidPCA Normalize --basename test_run_LOWMEM --path_to_directory ./
ResidPCA Standardize --basename test_run --path_to_directory ./
ResidPCA residPCA_fit --basename test_run_LOWMEM --path_to_directory ./
```

Python environment example command:
```
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
    lowmem=True
)

scExp.Normalize() 
scExp.Standardize()
```

Once you have instantiated your class and Normalized and Standardized you must **either** run ResidPCA or StandardPCA.

**For ResidPCA:**

Command line example  command:
```
ResidPCA residPCA_fit --basename test_run_LOWMEM --path_to_directory ./
```

Python environment example command:

```
scExp.residPCA_fit()
```

**StandardPCA:**

Command line example  command:
```
StandardPCA residPCA_fit --basename test_run_LOWMEM --path_to_directory ./
```

Python environment example command:


```
scExp.StandardPCA_fit()
```

** These are the same steps as the step by step directions, but in the lowmem setting, the package can only be run up until this point.


