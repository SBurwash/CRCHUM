import scprep



# Preprocess 10X

#' This function takes a data-frame (genes x cells), 
#' creates a Scprep object with it and filters the object for default tags such as  
#' Min and max nFeature_RNA and % of MT
#' 
#'
#' @param data data-frame
#' @param percent_mt integer [0-100]
#' @param max_features integer [0-Inf]
#' @param min_features integer [0-Inf]
#' @return Preprocessed Seurat object
#' @export

def preprocess_10X(data, name='10X-project', percent_mt=20, max_features=5000, min_features=200):
    
    #Remove empty cells and empty genes
    scprep_data = scprep.filter.filter_empty_cells(data)
    scprep_data = scprep.filter.filter_empty_genes(scprep_data)
    
    #Remove elements based on mythocondrial percentage
    mt_genes = scprep.select.get_gene_set(data, starts_with=["MT-", "mt-"])
    scprep_data = scprep.filter.filter_gene_set_expression(data=scprep_data, genes=mt_genes, percentile=100-percent_mt)
    
    #Remove elements based on number of cell and number of features
    scprep_data = scprep.filter.filter_library_size(data=scprep_data, cutoff=(min_features, max_features), keep_cells='between')
    scprep_data = scprep.filter.filter_rare_genes(data=scprep_data, cutoff=0, min_cells=3)

    return scprep_data
        


# Load 10X

# This function takes the path to a 10X output folder and instanciates the Scprep object
#
# @param file string (path to file)
# @param percent_mt integer [0-100]
# @param max_features integer [0-Inf]
# @param min_features integer [0-Inf]
# @return Preprocessed Seurat object
# @export

def load_10X(dir_path, name='10X-project',percent_mt=20, max_features=5000, min_features=200):
    data = scprep.io.load_10X(dir_path, sparse=True, gene_labels='both')
    return preprocess_10X(data, percent_mt, max_features, min_features)


# Reduce Dimensions

# This function reduces the dimensions of the Normalized Seurat object
# Runs PCA and then UMAP and then performs clustering

# @param seurat.obj S4 instance
# @param ndims integer [3-100]
# @param res double [0-3]
# @return Seurat object with reduction embeddings
# @export

def reduce_dim(data, ndims=15, res=.1):
    data_ln = scprep.normalize.library_size_normalize(data)
    data_sq = scprep.transform.sqrt(data_ln)

    return data_sq
