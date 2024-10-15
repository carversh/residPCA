
#!/usr/bin/env python
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

class condPCA(object):
    def __init__(self, 
                 count_matrix_path, 
                 object_columns, 
                 variable_genes_flavor = "seurat_v3",
                 metadata_path=None, 
                 vars_to_regress=True, 
                 n_PCs=200, 
                 random_seed=9989999, 
                 vargenes_IterPCA="all", 
                 vargenes_Stand_Cond="all", 
                 BIC=True, 
                 save_image_outputs = False, 
                 path_to_directory = "./", 
                 basename=f'CondPCA_run_{datetime.now().strftime("%Y-%m-%d_%H_%M_%S")}',  
                 global_ct_cutoff=0.2,
                 logged=False,
                 sparse_PCA=False):
        """
        Parameters
        ----------
        count_matrix:
            Count matrix that must be log-normalized and standardized

        metadata:
            metadata containing cell type labels named "celltype"

        object_columns:
            columns that will be one hot encoded/columns that are factors 

        vars_to_regress:
            list of variables to regress out

        """
        if save_image_outputs:
            
            # create directory to save files to, raise error if directory already exists
            directory_path = path_to_directory + basename # create directory path
            if not os.path.exists(directory_path):
                os.makedirs(directory_path)
                print(f"Directory '{directory_path}' created successfully.")
                
                self.directory_path = directory_path
            else:
                raise ValueError(f"Directory '{directory_path}' already exists.")
            
        self.count_matrix_path = count_matrix_path
        self.var_flavor = variable_genes_flavor
                
        count_matrix = scanpy.read(count_matrix_path) # cells x genes, pd.read_csv(count_matrix_path, sep='\t', header=0, index_col=0)
        
        self.vargenes_Stand_Cond = vargenes_Stand_Cond # number of variable genes for standard and conditional pca
        if self.vargenes_Stand_Cond == "all": # if all, make all genes are variable genes
            count_matrix.var['highly_variable'] = True
            print("added column")
        else:
            if self.var_flavor == "seurat_v3":
                try:
                    print("Finding most variable genes using 'seurat_v3' flavor.")
                    scanpy.pp.highly_variable_genes(count_matrix, flavor='seurat_v3', n_top_genes=vargenes_Stand_Cond)
                except Exception as e:
                    print(f"An error occurred: {e}")
                    print("Using 'seurat' flavor instead to compute variable genes on Standard and Conditional PCA.")
                    self.exception = True
            else:
                pass

        self.count_matrix = count_matrix.copy()

        # check if sparse counnt matrix and convert to dense if so
        if scipy.sparse.issparse(self.count_matrix.X):
            self.count_matrix.X = self.count_matrix.X.toarray()
     
        # read in h5ad of txt file
        if not vars_to_regress: # if vars_to_regress is False or an empty list
            self.metadata = False
        else:
            if count_matrix_path.endswith('.txt'):
                self.metadata = pd.read_csv(metadata_path, sep='\t', header=0, index_col=0)
            elif count_matrix_path.endswith('.h5ad'):
                self.metadata = self.count_matrix.obs.copy() 
            else:
                raise ValueError("Count matrix must be in .txt or .h5ad format.")
        
        if self.metadata is not False: # if there is metadata to regress out
            if isinstance(vars_to_regress, list):
                self.vars_to_regress = pd.Index(vars_to_regress)
            else:
                self.vars_to_regress = self.metadata.columns

        if not object_columns: # if no object columns are specified
            pass
        else:
            # one hot encode necessary metadata variables
            self.object_columns = object_columns # obtain columns that must be one hot encoded
            self.metadata[self.object_columns] = self.metadata[self.object_columns].astype(object) # convert these columns to objects
        
        self.random_seed = random_seed 
        np.random.seed(self.random_seed) # set random seed
        self.n_PCs = n_PCs # number of PCs to extract
        self.original_n_PCs = self.n_PCs # number of PCs to extract, self.n_PCs is modified in Iter_PCA_fit() function when not enough cells, and is edited back to the user designated number of PCs when necessary
        self.vargenes_IterPCA = vargenes_IterPCA # number of variable genes for iterative pca

        self.save_image_outputs = save_image_outputs # save image outputs
        self.BIC = BIC # compute BIC
        self.global_ct_cutoff  = global_ct_cutoff # cutoff of squared correlation in global vs ct specific states
        self.logged = logged # if the data is logged
        self.sparse_PCA = sparse_PCA # if performing sparse PCA instead of PCA

    def Normalize(self):
        """ 
        Normalize and take log1p of count data (for conditional and standard, not iterative)
        """
        # check if any rows or columns sum to 0 and throw error if so
        if np.any(self.count_matrix.X.sum(axis=1) == 0):
            raise ValueError("Error: Some Cells have zero counts.")
        if np.any(self.count_matrix.X.sum(axis=0) == 0):
            raise ValueError("Error: Some Genes have zero counts.")
        
        if self.logged:
            print("Data is already logged. Skipping log normalization.")
        else:
            # update scanpy object to normalize all rows, so every cell sums to 10k
            scanpy.pp.normalize_total(self.count_matrix, target_sum = 10000) 
            # log transform
            scanpy.pp.log1p(self.count_matrix) 

        # check whether alternative form of finding variable genes is necessary
        if hasattr(self, 'exception') and self.exception == True or self.var_flavor == "seurat":
            print("Finding most variable genes using 'seurat' flavor.")
            scanpy.pp.highly_variable_genes(self.count_matrix, n_top_genes=self.vargenes_Stand_Cond)
        self.count_matrix = self.count_matrix[:, self.count_matrix.var['highly_variable']] 

        # plot mean variance relationship if specified by user
        if self.save_image_outputs:
            self._plot_mean_variance_relationship(self.count_matrix.X, label="All Cells")
        

    def Standardize(self):
        """ 
        Standardize count data AND metadata (for conditional and standard, not iterative)
        """
        # Standardize count data
        # if some genes in counts matrix have zero standard deviation
        if np.any(np.std(self.count_matrix.X, axis=0) == 0):
            raise ValueError("Error: Some Genes have zero standard deviation.")
        # only subset the matrix to the most variable genes
        self.standardized_count_data = self._standardize(self.count_matrix.X)
        # Process metadata/covariates for standardization:
        if self.metadata is not False:
            # subset to only variables that you want to regress out
            self.metadata = self.metadata[self.vars_to_regress] 
            # WARNING IN FOLLOWING LINE BECAUSE CONVERTING OBJECT THAT LOOKS NUMERIC TO BE ONE HOT ENCODED, this is batch
            self.IterPCA_metadata = pd.get_dummies(self.metadata, drop_first=False)
            # Convert factor covariates to dummy variables dropping one column 
            self.metadata = pd.get_dummies(self.metadata, drop_first=True) 
            self.standardized_metadata = self._standardize(self.metadata)
    
    def _standardize(self, mat): # simple function performing standardization
        # compute means of genes or covariates
        mean_vector = np.mean(mat, axis=0)
       # compute standard deviation of genes or covariates
        std_vector = np.std(mat, axis=0)
        # standardize by gene or covariates
        stand_mat = (mat - mean_vector) / std_vector 
        return stand_mat
    
    def _regress_covariates(self, standardized_metadata, standardized_count_data): # function regressing set of covariates
        # append ones to standardized meta for intercept
        standardized_metadata = np.c_[np.ones((standardized_metadata.shape[0], 1)), standardized_metadata] 
        # compute inverse of np.matmul(A^T, A) where A is the standardized metadata or covariates
        #inv_cov = np.linalg.pinv(np.matmul(standardized_metadata.T, standardized_metadata) ) 
        inv_cov = np.linalg.inv(standardized_metadata.T @ standardized_metadata)
        # compute betas per gene
        betas = inv_cov @ standardized_metadata.T @ standardized_count_data
        # compute prediction
        #prediction = np.matmul(standardized_metadata, betas) # compute prediction
        prediction = standardized_metadata @ betas # compute prediction
        # compute residual
        residual = standardized_count_data - prediction 
        standardized_residual = self._standardize(residual)
        return standardized_residual
    
    def _fit_pca(self, mat, standardPCA, condPCA, iterPCA, iterPCA_genenames, iterPCA_cellnames, iterPCA_CT): # fitting PCA
        if self.sparse_PCA:
            # instantiate PCA with hyperparameters
            pca = SparsePCA(n_components=self.n_PCs, random_state=self.random_seed) 
        else:
            # instantiate PCA with hyperparameters
            pca = PCA(n_components=self.n_PCs, random_state=self.random_seed) 
        
        # projections (of input data onto eigenvectors)
        pca.fit(mat) 
        # retrieve eigenvectors/gene loadings
        gene_loadings = pca.components_ 
        # retrive cell embeddings
        cell_embeddings = pca.transform(mat)
        # retrieve eigenvalues
        eigenvalues = pca.explained_variance_        
        # if iterative PCA 
        if iterPCA:           
            # convert gene loadings to dataframe
            gene_loadings = pd.DataFrame(gene_loadings.T, index = list(iterPCA_genenames ), columns = [f'PC_{i}' for i in range(1, (gene_loadings.T.shape[1]+1))])               
            # convert cell embeddings to dataframe
            cell_embeddings = pd.DataFrame(cell_embeddings, index = list(iterPCA_cellnames), columns = [f'PC_{i}' for i in range(1, (cell_embeddings.shape[1]+1))])
        # convert eigenvalues to dataframe
        eigenvalues = pd.DataFrame(eigenvalues, index = [f'PC_{i}' for i in range(1, (eigenvalues.shape[0]+1))], columns=["Eigenvalues"])
        # if Standard or Conditional PCA, construct dataframes based on gene and cell list from original count matrix
        if not iterPCA: 
            # convert gene loadings to dataframe
            gene_loadings = pd.DataFrame(gene_loadings.T, index = list(self.count_matrix.var_names[self.count_matrix.var['highly_variable']] ), columns = [f'PC_{i}' for i in range(1, (gene_loadings.T.shape[1]+1))])               
            # convert cell embeddings to dataframe
            cell_embeddings = pd.DataFrame(cell_embeddings, index = list(self.count_matrix.obs_names), columns = [f'PC_{i}' for i in range(1, (cell_embeddings.shape[1]+1))])    
        if self.BIC:
            if standardPCA:
                # compute BIC
                min_BIC_index = self._compute_BIC(eigenvalues, self.standardized_residual, "Standard PCA")
                elbow_PC = self._compute_elbow(eigenvalues, self.standardized_residual, "Standard PCA")
            if condPCA:
                # compute BIC
                min_BIC_index = self._compute_BIC(eigenvalues, self.standardized_residual, "CondPCA")
                elbow_PC = self._compute_elbow(eigenvalues, self.standardized_residual, "CondPCA")
            if iterPCA:
                # compute BIC
                min_BIC_index = self._compute_BIC(eigenvalues, self.standardized_residual, "IterPCA", iterPCA_CT)
                elbow_PC = self._compute_elbow(eigenvalues, self.standardized_residual, "IterPCA", iterPCA_CT)
            BIC_cutoff = "PC_" + (min_BIC_index + 1).astype(str)
            return cell_embeddings, gene_loadings, eigenvalues, BIC_cutoff, elbow_PC
        return cell_embeddings, gene_loadings, eigenvalues, "Not Calculated", "Not Calculated"
    
    def _fit_model(self, standardized_metadata, standardized_count_data, standardPCA=False, condPCA = False, iterPCA=False, iterPCA_genenames=False, iterPCA_cellnames=False, iterPCA_CT=False): # regress out covariates and then input into PCA
        
        if standardized_metadata is not False: # if there is metadata to regress out
            # regress out covariates (including celltype) and retrieve standardized residual
            self.standardized_residual = self._regress_covariates(standardized_metadata = standardized_metadata, standardized_count_data= standardized_count_data) # REMOVE SELF  
        else: # no metadata, perform pca on standardized counts matrix
            print("No metadata to regress out.")
            self.standardized_residual = standardized_count_data # REMOVE SELF      
        # if not iterative PCA, able to add gene names and cell names here, but must subset if IterPCA
        if not iterPCA: 
            # return standardized residual as a dataframe with gene and cell names:
            self.standardized_residual = pd.DataFrame(self.standardized_residual, index = list(self.count_matrix.obs_names), columns = list(self.count_matrix.var_names[self.count_matrix.var['highly_variable']]))# REMOVE SELF
        if iterPCA:
            # return standardized residual as a dataframe with gene and cell names of the given subset:
            self.standardized_residual = pd.DataFrame(self.standardized_residual, index = list(iterPCA_cellnames), columns = list(iterPCA_genenames))# REMOVE SELF   
        
        # perform PCA on count matrix
        return self._fit_pca(self.standardized_residual, standardPCA=standardPCA, condPCA=condPCA, iterPCA=iterPCA, iterPCA_genenames=iterPCA_genenames, iterPCA_cellnames=iterPCA_cellnames, iterPCA_CT=iterPCA_CT)

    def _mapping_IterPCA_subset_dataframes_to_PCA(self, metadata, CT_exp_column): # function that subsets count matrix to a particular cell type and then performs PCA on that subset
        # remove "celltype_" from the string CT_exp_column
        CT_label = CT_exp_column.replace('celltype_', '')
        # extract indices of the cells that belong to the particular cell type of interest (indicated by CT_column, which is a column name)
        indices_given_ct = self.dataframe_CT_indices[CT_exp_column]
        # Check if the sum of indices is less than or equal to number of PCs requested
        if sum(indices_given_ct) <= self.n_PCs:
            print(f'Warning: Cell type {CT_label} has less than {sum(indices_given_ct)} cells and {self.n_PCs} principal components were requested.')       
            modified_num_PCs = int(sum(indices_given_ct) * 0.75)
            print(f'Reducing the number of PCs for {CT_label} Iterative PCA to {modified_num_PCs}.')
            # temporarily changing the number of PCs to extract just for this run and then resetting
            self.n_PCs = modified_num_PCs
        
        if metadata is not False: # if there is metadata to regress out
            # subset the count data to the cells belonging to the cell type
            metadata_subset_to_CT = metadata[indices_given_ct]
        # Re-process from log-normalized data to standadization of the matrix (find new set of variable genes for the subset)      
        # make a tmp copy and subset the matrix to cells in the particular cell type identify highly variable genes from log normalized count matrix
        # re-read in count data and subset
        count_matrix = scanpy.read(self.count_matrix_path) # cells x genes, pd.read_csv(count_matrix_path, sep='\t', header=0, index_col=0)
        count_matrix = count_matrix[indices_given_ct, :] # subset to cells in cell type

        #exception = False # marks whether variable genes method worked or not
        if self.vargenes_IterPCA == "all":
            count_matrix.var['highly_variable'] = True
        else:
            if self.var_flavor == "seurat_v3":
                try:
                    print(f'Using "seurat_v3" flavor to compute variable genes in Iterative PCA on {CT_label}.')
                    scanpy.pp.highly_variable_genes(count_matrix, flavor='seurat_v3', n_top_genes=self.vargenes_IterPCA)
                except Exception as e:
                    print(f"An error occurred: {e}")
                    print(f'Using "seurat" flavor to compute variable genes in Iterative PCA on {CT_label}.')
                    exception = True
            else:
                pass

        subset_iter = count_matrix.copy()

        if scipy.sparse.issparse(subset_iter.X):
            subset_iter.X = subset_iter.X.toarray()

        # update scanpy object to normalize all rows, so every cell sums to 10k
        scanpy.pp.normalize_total(subset_iter, target_sum = 10000)      
        # log transform
        scanpy.pp.log1p(subset_iter)  

        # if alternative way of finding variable genes is necessary because default method threw an error
        if 'exception' in locals() and exception == True or self.var_flavor == "seurat":
            print(f'Using "seurat" flavor to compute variable genes in Iterative PCA on {CT_label}.')
            scanpy.pp.highly_variable_genes(subset_iter, n_top_genes=self.vargenes_IterPCA)
        subset_iter = subset_iter[:, subset_iter.var['highly_variable']] 

        # plot mean variance relationship if specified by user
        if self.save_image_outputs:
            # if white space, replace with underscore
            CT_label = CT_label.replace(' ', '_')
            self._plot_mean_variance_relationship(subset_iter.X, label=CT_label)   

        # check if the counts matrix contains any genes with 0 variance
        if np.any(np.std(subset_iter.X, axis=0) == 0):
            if np.sum(np.std(subset_iter.X, axis=0) == 0) == subset_iter.X.shape[1]:
                # trigger an error if all genes have 0 variance
                raise KeyError("Error: All genes have 0 variance when performing iterative PCA on " + CT_label + " cells.")
            else:
                print("Warning: Some genes have 0 variance when performing iterative PCA on " + CT_label + " cells. " + str(np.sum(np.std(subset_iter.X, axis=0) == 0) ) + " genes will be removed from the " + CT_label + " counts matrix.")
                subset_iter = subset_iter[:, (~np.array(np.std(subset_iter.X, axis=0) == 0))]

        # Re-standardize count databecause it has just been subset
        log_norm_data_subset_to_CT = self._standardize(subset_iter.X)
        if metadata is not False: # if there is metadata to regress out  
            # Re-standardize metadata because it has just been subset
            metadata_subset_to_CT = self._standardize(metadata_subset_to_CT)
        # extract the cell names or barcodes of the cells belonging to the cell type of interest
        cellnames = subset_iter.obs_names 
        # extract the gene names of the genes belonging to the most variable genes within that subset
        genenames = subset_iter.var_names 
        if metadata is not False:
            # fit the given model by regressing out covariates and performing PCA
            output =  self._fit_model(standardized_metadata=metadata_subset_to_CT,standardized_count_data = log_norm_data_subset_to_CT, iterPCA=True, iterPCA_genenames=genenames, iterPCA_cellnames = cellnames, iterPCA_CT=CT_label)
        else:
            # fit the given model by regressing out covariates and performing PCA
            output =  self._fit_model(standardized_metadata=False,standardized_count_data = log_norm_data_subset_to_CT, iterPCA=True, iterPCA_genenames=genenames, iterPCA_cellnames = cellnames, iterPCA_CT=CT_label)
        # reset the number of PCs to the original input
        self.n_PCs = self.original_n_PCs
        
        
        checking = output + (pd.DataFrame(log_norm_data_subset_to_CT, index = list(cellnames), columns = list(genenames)),)# remove FOR CHECKING
        return output #checking## can remove final addition FOR CHECKING
    
    def CondPCA_fit(self): 
        if self.metadata is not False:     
            # fit linear model (regress out covariates) and fit PCA -- covariates contain cell type
            self.CondPCA_cell_embeddings, self.CondPCA_gene_loadings, self.CondPCA_eigenvalues, self.CondPCA_BIC_cutoff, self.CondPCA_elbow = self._fit_model(standardized_metadata=self.standardized_metadata,standardized_count_data= self.standardized_count_data, condPCA = True)
        else: 
            raise ValueError("Cannot perform Conditional PCA. No celltype column specified in metadata")

    def StandardPCA_fit(self):
        if self.metadata is not False: # if there is metadata to regress out
            if all(col.startswith('celltype_') for col in self.metadata): # if only metadata provided that is celltype
                print("Only celltype column provided in metadata. No covariates to regress out for Standard PCA.")
                self.StandardPCA_cell_embeddings, self.StandardPCA_gene_loadings, self.StandardPCA_eigenvalues, self.StandardPCA_BIC_cutoff, self.StandardPCA_elbow = self._fit_model(standardized_metadata=False,standardized_count_data= self.standardized_count_data,standardPCA=True)
            else:
                # remove celltype from covariate space
                standardized_metadata_minus_celltype = self.standardized_metadata.drop(columns = self.standardized_metadata.filter(like="celltype", axis=1).columns )
                # fit linear model (regress out covariates) and fit PCA -- covariates do not contain cell type
                self.StandardPCA_cell_embeddings, self.StandardPCA_gene_loadings, self.StandardPCA_eigenvalues, self.StandardPCA_BIC_cutoff, self.StandardPCA_elbow = self._fit_model(standardized_metadata=standardized_metadata_minus_celltype,standardized_count_data= self.standardized_count_data,standardPCA=True)
        else: # if there is no metadata to regress out
            self.StandardPCA_cell_embeddings, self.StandardPCA_gene_loadings, self.StandardPCA_eigenvalues, self.StandardPCA_BIC_cutoff = self._fit_model(standardized_metadata=False,standardized_count_data= self.standardized_count_data,standardPCA=True)


    def Iter_PCA_fit(self):
        if self.metadata is not False:   
            if all(col.startswith('celltype_') for col in self.metadata): # if only metadata provided that is celltype
                print("Only celltype column provided in metadata. No covariates to regress out for Iterative PCA.")
                pass

            else:
                # remove celltype from standardized covariate space
                metadata_minus_celltype = self.standardized_metadata.drop(columns = self.standardized_metadata.filter(like="celltype", axis=1).columns )  
            
            # get dataframe with boolean indices for each cell type
            self.dataframe_CT_indices = self.IterPCA_metadata.filter(like="celltype", axis=1).astype(bool)   
            # get the name of the columns that indicate a cell type
            celltype_colnames = self.dataframe_CT_indices.columns  
            # create empty dictionaries to store results of iterative PCA per cell type
            #self.IterPCA_residuals = {} # CAN REMOVE THIS IS USED FOR CHECKING
            self.IterPCA_cell_embeddings = {}
            self.IterPCA_gene_loadings = {}
            self.IterPCA_eigenvalues = {}
            self.IterPCA_BIC_cutoff = {}
            self.IterPCA_elbow = {}
            # iterate through each cell type and perform iterative PCA, storing results in dictionaries
            for celltype_column in celltype_colnames:
                # obtain cell type name, replace spaces with underscores
                tmp_CT = celltype_column.replace("celltype_", "").replace(" ", "_")
                if all(col.startswith('celltype_') for col in self.metadata): # if only metadata provided that is celltype
                    print("Only celltype column provided in metadata. No covariates to regress out for Iterative PCA.")
                    tmp_result = self._mapping_IterPCA_subset_dataframes_to_PCA(metadata=False, CT_exp_column=celltype_column)
                else: 
                    tmp_result = self._mapping_IterPCA_subset_dataframes_to_PCA(metadata_minus_celltype, celltype_column)
                # append results to appropriate dictionary
                self.IterPCA_cell_embeddings[tmp_CT] = tmp_result[0]
                self.IterPCA_gene_loadings[tmp_CT] = tmp_result[1]
                self.IterPCA_eigenvalues[tmp_CT] = tmp_result[2]
                self.IterPCA_BIC_cutoff[tmp_CT] = tmp_result[3]
                self.IterPCA_elbow[tmp_CT] = tmp_result[4]
                #self.IterPCA_residuals[tmp_CT]  = tmp_result[4] # CAN REMOVE USED FOR CHECKING
        else:
            raise ValueError("Cannot perform Iterative PCA. No celltype column specified in metadata")

    def _plot_mean_variance_relationship(self, log_normed_data, label):
        # compute the mean of every gene/column of the matrix
        mean_vector = np.mean(log_normed_data, axis=0)
        # compute the variance of every gene/column of the matrix
        var_vector = np.var(log_normed_data, axis=0)
        plt.scatter(mean_vector, var_vector)
        plt.xlabel("Mean")
        plt.ylabel("Variance")
        plt.title(f'Mean-Variance Relationship ({label}, {log_normed_data.shape[0]} Cells)')
        # save plot to current directory
        plt.savefig(f'{self.directory_path}/Mean_Variance_Relationship_{label}.png')
        # close plot
        plt.close()
    
    def _compute_elbow(self, eigenvalues, standardized_residual, label, ct=False):
        sum_eigenvalues = np.sum(eigenvalues)[0]
        # compute elbow
        prop_var_explained = np.zeros(eigenvalues.shape[0])
        cumulative_var_explained = np.zeros(eigenvalues.shape[0]) 
        for i in range(0, eigenvalues.shape[0]):
            prop_var_explained[i] = eigenvalues["Eigenvalues"][i] / sum_eigenvalues
            cumulative_var_explained[i] = np.sum(prop_var_explained[0:i+1]) 

        # find where the cumulative variance explained begins decreasing by less than 0.05
        indices = np.where(np.abs(np.diff(prop_var_explained)) < 0.001)[0] # REMOVE SELF

        # Find the first index
        first_index = indices[0] if len(indices) > 0 else None

        elbow = eigenvalues.iloc[first_index].name
        
        if self.save_image_outputs:
            # Plot BIC values
            plt.plot(range(1, len(eigenvalues["Eigenvalues"][0:50]) + 1), eigenvalues["Eigenvalues"][0:50], marker='o', linestyle='-')
            plt.xlabel('Principal Component Number')
            plt.ylabel('Eigenvalue')
            plt.axvline(x=first_index + 1, color='r', linestyle='--')
            plt.legend(['Eigenvalues', elbow])
            plt.grid(True)
            if not ct:
                plt.title(f'{label} Eigenvalues vs. Number of Principal Components')
                plt.savefig(f'{self.directory_path}/Elbow_plot_{label.replace(" ", "_")}.png')
            if ct:
                plt.title(f'{label} Eigenvalues vs. Number of Principal Components ({ct})')
                plt.savefig(f'{self.directory_path}/Elbow_plot_{label.replace(" ", "_")}_{ct.replace(" ", "_")}.png')
            plt.close()
        return elbow

    def _compute_BIC(self, eigenvalues, standardized_residual, label, ct=False):
        # extract the standardized residuals as a numpy array
        X = standardized_residual.values
        # Compute the sample covariance matrix
        cov_matrix = (X.T @ X) / (X.shape[0] - 1)
        # Compute the trace of the sample covariance matrix (this is equal to the sum of the eigenvalues)
        trace = np.trace(cov_matrix)
        # compute BIC
        p = X.shape[1] 
        n = X.shape[0] 
        # Initialize an array to store BIC values
        PCs = eigenvalues.shape[0]
        BIC_values = np.zeros(PCs)
        # Perform calculations for each value of j
        for j in range(0, PCs ):
            ell_bar = (trace - np.sum(eigenvalues.iloc[0:j, 0]) ) / (p - j)
            epsilon = 1e-10
            ell_bar = epsilon if ell_bar <= 0 else ell_bar # set ell_bar to a value very close to 0 if it's negative
            dj = (j + 1) * (p + 1 - j / 2)
            term_1 = n * np.log(np.prod(eigenvalues.iloc[0:j, 0]))
            term_2 = n * (p - j) * np.log(ell_bar)
            term_3 = np.log(n) * dj
            term_4 = n * (np.log((n - 1) / n )**p) + n * p * (1 + np.log(2 * np.pi))
            BIC_j = term_1 + term_2 + term_3 + term_4
            # Store BIC value in the array
            BIC_values[j] = BIC_j
        # Find the index corresponding to the minimum BIC value
        min_BIC_index = np.argmin(BIC_values)
        if self.save_image_outputs:
            # Plot BIC values
            plt.plot(range(1, PCs + 1), BIC_values, marker='o', linestyle='-', color='b')
            plt.axvline(x=min_BIC_index + 1, color='r', linestyle='--', label=f'Min BIC at $j = {min_BIC_index + 1}$')
            plt.xlabel('Number of Principal Components (j)')
            plt.ylabel('BIC Value')           
            plt.legend()
            plt.grid(True)
            if not ct:
                plt.title(f'{label} BIC Values vs. Number of Principal Components')
                plt.savefig(f'{self.directory_path}/BIC_plot_{label.replace(" ", "_")}.png')
            if ct:
                plt.title(f'{label} BIC Values vs. Number of Principal Components ({ct})')
                plt.savefig(f'{self.directory_path}/BIC_plot_{label.replace(" ", "_")}_{ct.replace(" ", "_")}.png')
            plt.close()
        return min_BIC_index
    
    # function that subsets a dataframe to the column that is the BIC cutoff
    def _sub_dataframe_BIC(self, gene_loadings, BIC_cuttoff):
        stop_variable_index = gene_loadings.columns.get_loc(BIC_cuttoff) + 1
        gene_loadings_sub_BIC = gene_loadings.iloc[:, :stop_variable_index]
        return gene_loadings_sub_BIC

    # function that computes the correlation of each column of dataframe 1 with each column of dataframe 2
    def _compute_correlation(self, df1, df2):
        # create empty dataframe to store results
        df = pd.DataFrame(index=df1.columns, columns=df2.columns)
        # iterate through each column of dataframe 1
        for col1 in df1.columns:
            # iterate through each column of dataframe 2
            for col2 in df2.columns:
                # compute correlation between the two columns
                df.loc[col1, col2] = np.corrcoef(df1[col1], df2[col2])[0,1]
        return df

    # function that accepts a Standard/CondPCA gene loadings dataframe and dictionary of IterPCA gene loadings dataframes and returns the squared correlation between the Standard/CondPCA gene loadings dataframe and each of the IterPCA gene loadings dataframes in the dictionary
    def _compute_squared_correlation_w_iterative(self, gene_loadings, gene_loadings_dict):       
        # initialize empty dataframe to store results
        squared_correlations = pd.DataFrame(index=gene_loadings.columns, columns=gene_loadings_dict.keys())
        # iterate through each cell type and compute the squared correlation between the gene loadings and the gene loadings for that cell type
        for celltype in gene_loadings_dict.keys():
            # computing intersecting genes between two dataframes and subset both dataframes by those genes
            intersecting_gene_names = gene_loadings.index.intersection(gene_loadings_dict[celltype].index)
            gene_loadings_sub = gene_loadings.loc[intersecting_gene_names,]
            gene_loadings_dict_sub = gene_loadings_dict[celltype].loc[intersecting_gene_names,]
            # compute correlation between the two dataframes and square it
            corr = self._compute_correlation(gene_loadings_sub, gene_loadings_dict_sub)**2
            # find the maximum correlation for each row, this will compute the maximum correlation in each PC of Standard of CondPCA
            max_corr = corr.max(axis=1)
            # append squared correlation to dataframe
            squared_correlations[celltype] = max_corr
        # return the list of squared correlations
        return squared_correlations
    
    # function that accepts Standard/CondPCA gene loadings dataframe and depending on order, outputs a dataframe of the cross correlation between those two datasets, the first argument will retain its dimensionality in number of PCs
    def _compute_squared_correlation_w_StandCond(self, gene_loadings1, gene_loadings2, gene_loadings2_label):
        # initialize empty dataframe to store results
        squared_correlations = pd.DataFrame(index=gene_loadings1.columns, columns=[gene_loadings2_label])
        # computing intersecting genes between two dataframes and subset both dataframes by those genes
        intersecting_gene_names = gene_loadings1.index.intersection(gene_loadings2.index)
        gene_loadings1_sub = gene_loadings1.loc[intersecting_gene_names,]
        gene_loadings2_sub = gene_loadings2.loc[intersecting_gene_names,]
        # compute correlation between the two dataframes and square it
        corr = self._compute_correlation(gene_loadings1_sub, gene_loadings2_sub)**2
        # find the maximum correlation for each row, this will compute the maximum correlation in each PC of Standard of CondPCA
        squared_correlations[gene_loadings2_label] = corr.max(axis=1)
        return squared_correlations
    
    def _plot_heatmap_global_ct_specific(self, data, PCA_type):
        # Create a mask for values greater than designated
        mask = data > self.global_ct_cutoff
        plt.figure(figsize=(12, 13))
        sns.heatmap(data, annot=True, cmap='Reds', linewidths=.5, vmin=0, vmax=1, cbar=False, fmt=".2f")
        # Add boxes only around entries greater than 0.20
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                if mask.iloc[i, j]:
                    rect = plt.Rectangle((j, i), 1, 1, linewidth=1, edgecolor='black', facecolor='none')
                    plt.gca().add_patch(rect)
        plt.title(f'Heatmap of Squared Correlations Between {PCA_type} and Iterative PCA\n(only values greater than {self.global_ct_cutoff} are boxed)')
        plt.savefig(f'{self.directory_path}/heatmap_squared_correlations_betw_{PCA_type.replace(" ", "")}_and_IterPCA.png')
        plt.close()

    def _label_Global_CellType_States(self, df_squared_correlations, standard_or_conditional):
        # Goal: add CT_involved column, which tells you which cell types are involved in a state and Global_vs_CT column, which tells you if the state global or cell type specific
        # obtain boolean of where squared correlation is greater than cutoff
        greater_than_cutoff = df_squared_correlations > self.global_ct_cutoff
        # sum across rows to see if any entries are greater than cutoff or all are less than cutoff
        sum_rows_greater_than_cutoff = greater_than_cutoff.sum(axis=1)

        sum_rows_greater_than_cutoff = sum_rows_greater_than_cutoff.copy()

        # Create a new array of labels based on conditions
        state_labels = np.where(sum_rows_greater_than_cutoff == 0, f'{standard_or_conditional}-Only State',
                        np.where(sum_rows_greater_than_cutoff == greater_than_cutoff.shape[1], "All Cell Types State", "Cell Type Specific State"))
        # Function to find columns with True values in each row
        def true_columns(row):
            return greater_than_cutoff.columns[row].tolist()
        # Apply the function to each row
        greater_than_cutoff['True_Columns'] = greater_than_cutoff.apply(lambda row: true_columns(row), axis=1)
        output = pd.DataFrame({"CT_involved": greater_than_cutoff['True_Columns']})
        output['Global_vs_CT'] = state_labels
        return output
    
    # Calculate proportions global vs ct specific
    def _calc_prop(self,df, conditional_or_standard):
        proportions = df['Global_vs_CT'].value_counts(normalize=True).to_frame(name='Proportion')
        proportions.columns =  [conditional_or_standard + ' PCA']
        return proportions
    
    # Plot proportions global vs ct specific
    def _prop_plot(self,final_df):
        # Plotting a stacked bar plot
        ax = final_df.T.plot(kind='bar', stacked=True)
                
        # Add Title and Labels
        plt.title('State Type Distribution')
        # Remove legend title
        ax.legend(title='')
        plt.xlabel('')
        plt.xticks(rotation=0)
        ax.legend(title='', bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.ylabel('Proportion')
        plt.savefig(f'{self.directory_path}/barplot_cell_state_distributions.png', bbox_inches='tight')
        plt.close()

    # create Upset Plot of cell type specific states
    def _Upset(self,df,label):    
        # convert data to be formatted for Upset plot
        data = df['CT_involved'].to_list()

        # Create UpSet plot
        example =from_memberships(data, data=np.arange((len(data)*7)).reshape(len(data), 7))
        plot(example, show_counts=True)
        plt.suptitle(f'{label} PCA - Upset Plot of Cell Type Specific States')
        plt.savefig(f'{self.directory_path}/{label}_Upset_plot_states.png', bbox_inches='tight')
        plt.close()
    
    # plot distributions of max squared correlations across all cell types
    def _plot_KDE(self, data, PCA_type):
        sns.kdeplot(data, fill=True, clip=(0, None))
        plt.xlabel('Squared Correlation')
        plt.title(f'KDE Plot of Max Corr.^2 between {PCA_type} and Iterative PCA, BIC Significant States')
        label=PCA_type.replace(" ", "_")
        plt.savefig(f'{self.directory_path}/{label}_cell_type_distribution.png', bbox_inches='tight')
        plt.close()
    
    def ID_Global_CellType_States(self):
        # subset to BIC cutoff if dataset is present
        if hasattr(self, 'StandardPCA_gene_loadings'):
            standard_gene_loadings_sub_BIC = self._sub_dataframe_BIC(self.StandardPCA_gene_loadings, self.StandardPCA_BIC_cutoff)
        if hasattr(self, 'CondPCA_gene_loadings'):
            cond_gene_loadings_sub_BIC = self._sub_dataframe_BIC(self.CondPCA_gene_loadings, self.CondPCA_BIC_cutoff)
        if hasattr(self, 'IterPCA_gene_loadings'):
                iter_gene_loadings_sub_BIC = {}
                for celltype in self.IterPCA_gene_loadings.keys():
                    iter_gene_loadings_sub_BIC[celltype] = self._sub_dataframe_BIC(self.IterPCA_gene_loadings[celltype], self.IterPCA_BIC_cutoff[celltype])  

        # compute cross correlations
        if hasattr(self, 'StandardPCA_gene_loadings') and hasattr(self, 'CondPCA_gene_loadings') and hasattr(self, 'IterPCA_gene_loadings'):
            # if Standard, Conditional, and Iterative exist
            
            # compute the squared correlation between the gene loadings for standard PCA and the gene loadings for each cell type
            StandardPCA_IterPCA_squared_correlations = self._compute_squared_correlation_w_iterative(standard_gene_loadings_sub_BIC, iter_gene_loadings_sub_BIC).apply(pd.to_numeric, errors='coerce') # RENAME MORAB

            # compute the squared correlation between the gene loadings for standard PCA and the gene conditional PCA
            StandardPCA_CondPCA_squared_correlations = self._compute_squared_correlation_w_StandCond(standard_gene_loadings_sub_BIC, cond_gene_loadings_sub_BIC, "CondPCA").apply(pd.to_numeric, errors='coerce')

            # combine standard correlations into one dataframe
            StandardPCA_correlations = pd.concat([StandardPCA_IterPCA_squared_correlations, StandardPCA_CondPCA_squared_correlations], axis=1)


            # compute the squared correlation between the gene loadings for standard PCA and the gene loadings for each cell type
            CondPCA_IterPCA_squared_correlations = self._compute_squared_correlation_w_iterative(cond_gene_loadings_sub_BIC, iter_gene_loadings_sub_BIC).apply(pd.to_numeric, errors='coerce')

            # compute the squared correlation between the gene loadings for conditional PCA and the gene loadings for standard PCA
            CondPCA_StandardPCA_squared_correlations = self._compute_squared_correlation_w_StandCond(cond_gene_loadings_sub_BIC, standard_gene_loadings_sub_BIC, "StandardPCA").apply(pd.to_numeric, errors='coerce')

            # combine conditional correlations into one dataframe
            CondPCA_correlations = pd.concat([CondPCA_IterPCA_squared_correlations, CondPCA_StandardPCA_squared_correlations], axis=1)

            # label states with their corresponding cell type
            self.StandardPCA_correlations = pd.merge(StandardPCA_correlations, self._label_Global_CellType_States(StandardPCA_correlations.drop(["CondPCA", ], axis=1), "Standard"), left_index=True, right_index=True)
            self.CondPCA_correlations = pd.merge(CondPCA_correlations, self._label_Global_CellType_States(CondPCA_correlations.drop(["StandardPCA", ], axis=1), "Conditional"), left_index=True, right_index=True) 

            # COMPUTE IMAGE OUTPUTS
            if self.save_image_outputs:
                # plot heatmaps
                self._plot_heatmap_global_ct_specific(StandardPCA_correlations, "StandardPCA")
                self._plot_heatmap_global_ct_specific(CondPCA_correlations, "CondPCA")
                # Calculate proportions global vs ct specific states
                standard_prop = self._calc_prop(self.StandardPCA_correlations.drop(["CondPCA", ], axis=1), "Standard")
                conditional_prop = self._calc_prop(self.CondPCA_correlations.drop(["StandardPCA", ], axis=1), "Conditional")
                # Outer combine DataFrames and impute 0 for missing values
                combined_prop = pd.merge(standard_prop, conditional_prop, left_index=True, right_index=True, how='outer').fillna(0)
                self._prop_plot(combined_prop)
                # also make Upset plots
                self._Upset(self.StandardPCA_correlations.drop(["CondPCA", ], axis=1), "Standard")
                self._Upset(self.CondPCA_correlations.drop(["StandardPCA", ], axis=1), "Conditional")
                self._plot_KDE(StandardPCA_IterPCA_squared_correlations, PCA_type="Standard PCA")
                self._plot_KDE(CondPCA_IterPCA_squared_correlations, PCA_type="Conditional PCA")        
                

        elif hasattr(self, 'StandardPCA_gene_loadings') and hasattr(self, 'CondPCA_gene_loadings'):
            # if only standard and conditional exist
            
            # compute the squared correlation between the gene loadings for standard PCA and the gene conditional PCA
            self.StandardPCA_correlations = self._compute_squared_correlation_w_StandCond(standard_gene_loadings_sub_BIC, cond_gene_loadings_sub_BIC, "CondPCA").apply(pd.to_numeric, errors='coerce')
            
            # compute the squared correlation between the gene loadings for conditional PCA and the gene loadings for standard PCA
            self.CondPCA_correlations = self._compute_squared_correlation_w_StandCond(cond_gene_loadings_sub_BIC, standard_gene_loadings_sub_BIC, "StandardPCA").apply(pd.to_numeric, errors='coerce')
            
        elif hasattr(self, 'StandardPCA_gene_loadings') and hasattr(self, 'IterPCA_gene_loadings'):
            # if only standard and iterative exist

            # compute the squared correlation between the gene loadings for standard PCA and the gene loadings for each cell type
            StandardPCA_correlations = self._compute_squared_correlation_w_iterative(standard_gene_loadings_sub_BIC, iter_gene_loadings_sub_BIC).apply(pd.to_numeric, errors='coerce')  

            # label states with their corresponding cell type
            self.StandardPCA_correlations = pd.merge(StandardPCA_correlations, self._label_Global_CellType_States(StandardPCA_correlations, "Standard"), left_index=True, right_index=True) # CHANGE FUNCTION
            
            # COMPUTE IMAGE OUTPUTS
            if self.save_image_outputs:
                # plot heatmaps
                self._plot_heatmap_global_ct_specific(StandardPCA_correlations, "StandardPCA")
                # Calculate proportions global vs ct specific states
                standard_prop = self._calc_prop(self.StandardPCA_correlations, "Standard")
                self._prop_plot(standard_prop)
                # also make Upset plots
                self._Upset(self.StandardPCA_correlations, "Standard")
                self._plot_KDE(StandardPCA_IterPCA_squared_correlations, PCA_type="Standard PCA")

        elif hasattr(self, 'CondPCA_gene_loadings') and hasattr(self, 'IterPCA_gene_loadings'):
            #if only conditional and iterative exist
            
            # compute the squared correlation between the gene loadings for standard PCA and the gene loadings for each cell type
            CondPCA_correlations = self._compute_squared_correlation_w_iterative(cond_gene_loadings_sub_BIC, iter_gene_loadings_sub_BIC).apply(pd.to_numeric, errors='coerce')

            # label states with their corresponding cell type
            self.CondPCA_correlations = pd.merge(CondPCA_correlations, self._label_Global_CellType_States(CondPCA_correlations, "Conditional"), left_index=True, right_index=True) # CHANGE FUNCTION

            # COMPUTE IMAGE OUTPUTS
            if self.save_image_outputs:
                # plot heatmaps
                self._plot_heatmap_global_ct_specific(CondPCA_correlations, "CondPCA")
                # Calculate proportions global vs ct specific states
                conditional_prop = self._calc_prop(self.CondPCA_correlations, "Conditional")
                self._prop_plot(conditional_prop)
                # also make Upset plots
                self._Upset(self.CondPCA_correlations, "Conditional")
                self._plot_KDE(CondPCA_IterPCA_squared_correlations, PCA_type="Conditional PCA")  

        elif hasattr(self, 'StandardPCA_gene_loadings'):
            raise ValueError("Only Standard PCA has been performed, must perform both/either Conditional or Iterative PCA as well.")
        elif hasattr(self, 'CondPCA_gene_loadings'):
            raise ValueError("Only Conditional PCA has been performed, must perform both/either Standard or Iterative PCA as well.")
        elif hasattr(self, 'IterPCA_gene_loadings'):
            raise ValueError("Only Iterative PCA has been performed, must perform both/either Standard or Conditional PCA as well.")
        else:
            raise ValueError("No processed datasets to compare.")

def main():
    pass   

if __name__=="__main__":
    main()