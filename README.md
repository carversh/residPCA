# residPCA Package v1.0.0

condPCA is an analysis pipeline for inferring cell states within and across cell types by leveraging cell type labels. 

It takes a count matrix (N cells x G genes) as input and produces ...

![Inline Image](https://github.com/carversh/residPCA/blob/main/pipeline_visual.png)

# Installation

```
pip install residPCA
```

# Step by Step Guide 

### Step 1 - instantiate a class with your input data

Example command:

```
scExp = condPCA(count_matrix_path="matrix.txt", metadata_path="metadata.txt", object_columns=['Batch', 'Sex','celltype'], save_image_outputs=False, BIC=True)

```

Input Data
  - ```count_matrix_path``` - points to the tab delimited file that is cells by genes dimensional, where rownames correspond to the barcoded cells and columnnames correspond to the genes, or a Scanpy object.
  - ```metadata_path``` - points to the tab delimited file that is cells by covariates/features dimensional.There must be a single column with a column name "celltype" that contains the cell type labels corresponding to each barcoded cell in the count matrix.

Parameters
  - ```object_columns = []``` - A list that specifies column names whose column should be treated as a factor or object and not a numeric column. The method will one hot encode these columns. Celltype must be in this list. Commonly, batch is another covariate that is included in this list.
  - ```n_PCs``` - number of PCs to output per method. Default is 200.
  - ```random_seed``` - Default is 998999.
  - ```vargenes_Stand_Cond``` - The number of variable genes to analyze during Standard and Conditional PCA. This can be set to an integer or "all" if all genes should be considered for analyses.
  - ```vargenes_IterPCA``` - The number of variable genes to analyze during Iterative PCA. This can be set to an integer or "all" if all genes should be considered for analyses.
  - ```BIC``` - True/False. Specified whether to compute the Bayesian Information Criterion after each method so that the set of states has a statistical cutoff. 
  - ```save_image_outputs``` - True/False. Specifies whether to save the graphical outputs during the pipeline processing (i.e. Mean-Variance Log-Normalization Graphs, BIC Cutoff Graphs,...).

  - ```path_to_directory``` - path where the basename directory will be created that will contain all images and data created.
  - ```basename``` - basename of the directory created that will contain all images and data created from each command. Default is ```CondPCA_run_{date}```, i.e. ```CondPCA_run_2024-01-05_13_35_14```.
  - ```global_ct_cutoff``` - squared correlation at which to above this cutoff, a cell type is involved in this given state and below this curoff, a cell type is not involved in this given state. Default is 0.2. (used in ```.ID_Global_CellType_States()``` function)
### Step 2 - log-normalize the count data

Example command:

```
scExp.Normalize()
```

### Step 3 - standardize the count data

Example command:

```
scExp.Standardize()
```

### Step 4 - perform Standard PCA

Example command:

```
scExp.StandardPCA_fit()
```
Returns the Standard PCA output in the form of dataframes. 

Outputs
  - ```scExp.StandardPCA_cell_embeddings``` - cell embeddings outputted by Standard PCA
  - ```scExp.StandardPCA_gene_loadings``` - gene loadings or eigenvectors outputted by Standard PCA
  - ```scExp.StandardPCA_eigenvalues``` - eigenvalues outputted by Standard PCA
  - ```scExp.StandardPCA_BIC_cutoff``` - PC cutoff that specifies the maximum state that is significant. For significant states, subset the cell embeddings and gene loadings from PC1 to the PC specified in this variable

### Step 5 - perform Conditional PCA (CondPCA)

Example command:

```
scExp.CondPCA_fit()
```
Returns the CondPCA output in the form of dataframes. 

Outputs
  - ```scExp.CondPCA_cell_embeddings``` - cell embeddings outputted by Conditional PCA
  - ```scExp.CondPCA_gene_loadings``` - gene loadings or eigenvectors outputted by Conditional PCA
  - ```scExp.CondPCA_eigenvalues``` - eigenvalues outputted by Conditional PCA
  - ```scExp.CondPCA_BIC_cutoff``` - PC cutoff that specifies the maximum state that is significant. For significant states, subset the cell embeddings and gene loadings from PC1 to the PC specified in this variable

### Step 6 - perform Iterative PCA (IterPCA)

Example command:

```
scExp.Iter_PCA_fit()
```
Returns the IterPCA output in the form of dictionaries, where each dictionary is equal to the length of the number of cell types. The keys of the dictionary correspond to the cell type in the "celltype" column of the metadata while the values of the dictionary represent the respective dataframe that corresponds to that cell type. 

Outputs
  - ```scExp.IterPCA_cell_embeddings``` - dictionary containing the cell embeddings outputted by Iterative PCA per cell type
  - ```scExp.IterPCA_gene_loadings``` - dictionary containing the gene loadings outputted by Iterative PCA per cell type
  - ```scExp.IterPCA_eigenvalues``` - dictionary containing the gene eigenvalues outputted by Iterative PCA per cell type
  - ```scExp.IterPCA_BIC_cutoff``` - dictionary of PC cutoffs that specifies the maximum state that is significant per cell type. For significant states, subset the cell embeddings and gene loadings from PC1 to the PC specified in this variable

Warning
  - If there are fewer than 200 cells in a cell type, the method will return an empty dataset for that given cell type.

### Step 7 - identify states that are cell type specific and states that span all cell types

Example command:

```
scExp.ID_Global_CellType_States()
```

Returns a dataframe for StandardPCA based states and CondPCA based states. Returns the maximum squared correlation between states identified in IterPCA and StandardPCA/CondPCA. Annotates each state and whether it is a global or cell type specific state as well as the cell types that are in that state.

Note: This method labels each state identified in StandardPCA and CondPCA as a global state, or a state spanning multiple cell types, or cell type specific states, or a state within a specific cell type. Additionally, this method labels which states belong to which cell types. To perform this method, ```scExp.Iter_PCA_fit()``` must be run beforehand and at least ```scExp.StandardPCA_fit()``` or ```scExp.CondPCA_fit()```  must be run. If only ```scExp.CondPCA_fit()``` is run, only the states belonoging to CondPCA will be evaluated, if only ```scExp.StandardPCA_fit()``` is run, only the states belonging to StandardPCA will be evaluated, if both are run, both will be evaluated.

Outputs 
 -```scExp.StandardPCA_IterPCA_squared_correlations``` - dataframe containing each state from StandardPCA and a row and the maximum correlation with a certain run of IterPCA (columns are labeled by IterPCA run on that given cell type). Two additional columns are included called "CT_involved", which tells you which cell types are involved in the state, and "Global_vs_CT", which tell you whether the state is global or cell type specific.
 - ```scExp.CondPCA_IterPCA_squared_correlations``` - dataframe containing each state from CondPCA and a row and the maximum correlation with a certain run of IterPCA (columns are labeled by IterPCA run on that given cell type). Two additional columns are included called "CT_involved", which tells you which cell types are involved in the state, and "Global_vs_CT", which tell you whether the state is global or cell type specific.

 If ```save_image_outputs == True```:
   - example

# Image Outputs

When initializing your experiment, if you set ```save_image_outputs = True```, a directory will be created that will contain all image based outputs from the method. Different commands will yield different image outputs depending on the task of the command. All relevant images and data will be saved to a directory called ```basename``` with the initial appended path ```path_to_directory```.
