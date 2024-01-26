#!/usr/bin/env Rscript

# Author_and_contribution: Niklas Mueller-Boetticher; created template
# Author_and_contribution: Lucie Pfeiferova; functions for Giotto HMRF spatial domain exploring
# Author_and_contribution: Søren Helweg Dam; created environment setup script, updated environment yaml, added configs and code for remaining cluster functions, tidied code

suppressPackageStartupMessages({
    library(optparse)
    library(jsonlite)
    library(SpatialExperiment)
    library(Giotto)
})

option_list <- list(
  make_option(
    c("-c", "--coordinates"),
    type = "character", default = NULL,
    help = "Path to coordinates (as tsv)."
  ),
  make_option(
    c("-m", "--matrix"),
    type = "character", default = NA,
    help = "Path to (transformed) counts (as mtx)."
  ),
  make_option(
    c("-f", "--features"),
    type = "character", default = NULL,
    help = "Path to features (as tsv)."
  ),
  make_option(
    c("-o", "--observations"),
    type = "character", default = NULL,
    help = "Path to observations (as tsv)."
  ),
  make_option(
    c("-n", "--neighbors"),
    type = "character", default = NA,
    help = "Path to neighbor definitions. Square matrix (not necessarily symmetric) where each row contains the neighbors of this observation (as mtx)."
  ),
  make_option(
    c("-d", "--out_dir"),
    type = "character", default = NULL,
    help = "Output directory."
  ),
  make_option(
    c("--dim_red"),
    type = "character", default = NA,
    help = "Reduced dimensionality representation (e.g. PCA)."
  ),
  make_option(
    c("--image"),
    type = "character", default = NA,
    help = "Path to H&E staining."
  ),
  make_option(
    c("--n_clusters"),
    type = "integer", default = NULL,
    help = "Number of clusters to return."
  ),
  make_option(
    c("--technology"),
    type = "character", default = NULL,
    help = "The technology of the dataset (Visium, ST, imaging-based)."
  ),
  make_option(
    c("--seed"),
    type = "integer", default = NULL,
    help = "Seed to use for random operations."
  ),
  make_option(
    c("--config"),
    type = "character", default = NA,
    help = "Optional config file (json) used to pass additional parameters."
  )
)
## -- make option for betas (eg c(0,8,10) and betas_to_add (eg which one to use)
description <- "Giotto: Implementations of HMRF, Leiden, Louvain, randomwalk, SNNclust, kmeans, and hierarchical"

opt_parser <- OptionParser(
  usage = description,
  option_list = option_list
)
opt <- parse_args(opt_parser)

out_dir <- opt$out_dir

# Output files
label_file <- file.path(out_dir, "domains.tsv")
embedding_file <- file.path(out_dir, "embedding.tsv")
# if additional output files are required write it also to out_dir

# Use these filepaths as input ...
coord_file <- opt$coordinates
feature_file <- opt$features
observation_file <- opt$observations

if (!is.na(opt$neighbors)) {
  neighbors_file <- opt$neighbors
}
if (!is.na(opt$matrix)) {
  matrix_file <- opt$matrix
}
if (!is.na(opt$dim_red)) {
  dimred_file <- opt$dim_red
}
if (!is.na(opt$image)) {
  image_file <- opt$image
}
if (!is.na(opt$config)) {
  config_file <- opt$config
  config <- fromJSON(config_file)
}

technology <- opt$technology
n_clusters <- opt$n_clusters

# You can get SpatialExperiment directly
get_SpatialExperiment <- function(
    feature_file,
    observation_file,
    coord_file,
    matrix_file = NA,
    reducedDim_file = NA,
    assay_name = "counts",
    reducedDim_name = "reducedDim") {
  rowData <- read.delim(feature_file, stringsAsFactors = FALSE, row.names = 1)
  colData <- read.delim(observation_file, stringsAsFactors = FALSE, row.names = 1)

  coordinates <- read.delim(coord_file, sep = "\t", row.names = 1)
  coordinates <- as.matrix(coordinates[rownames(colData), ])
  coordinates[,c(1:2)] <- as.numeric(coordinates[,c(1:2)])
    
  spe <- SpatialExperiment::SpatialExperiment(
    rowData = rowData, colData = colData, spatialCoords = coordinates
  )

  if (!is.na(matrix_file)) {
    assay(spe, assay_name, withDimnames = FALSE) <- as(Matrix::t(Matrix::readMM(matrix_file)), "CsparseMatrix")
  }

  # Filter features and samples
  if ("selected" %in% colnames(rowData(spe))) {
    spe <- spe[as.logical(rowData(spe)$selected), ]
  }
  if ("selected" %in% colnames(colData(spe))) {
    spe <- spe[, as.logical(colData(spe)$selected)]
  }

  if (!is.na(reducedDim_file)) {
    dimRed <- read.delim(reducedDim_file, stringsAsFactors = FALSE, row.names = 1)
    reducedDim(spe, reducedDim_name) <- as.matrix(dimRed[colnames(spe), ])
  }
  return(spe)
}


# Seed
seed <- opt$seed
set.seed(seed)

# SpatialExperiment
spe <- get_SpatialExperiment(
    feature_file = feature_file,
    observation_file = observation_file,
    coord_file = coord_file,
    matrix_file = matrix_file,
    reducedDim_file = dimred_file
)

## Configuration
method <- config$method
k <- config$k
type <- config$type

## Giotto instructions
python_path <- Sys.which(c("python"))
instrs <- createGiottoInstructions(save_dir = out_dir,
                                  save_plot = TRUE,
                                  show_plot = TRUE,
                                  python_path = python_path)



## raw expression counts expected
createGiotto_fn = function(spe, annotation = FALSE, selected_clustering = NULL, instructions = NULL){
  raw_expr <- SummarizedExperiment::assay(spe, "counts")
  #colnames(raw_expr) <- colData(sce)[,"Barcode"]
  #norm_expression <- SummarizedExperiment::assay(sce, "logcounts")
  
  cell_metadata <- SingleCellExperiment::colData(spe)
  cell_metadata$cell_ID <- rownames(SingleCellExperiment::colData(spe))
  colnames(cell_metadata)[c(1,2)] <- c("sdimx", "sdimy")
  cell_metadata <- as.data.frame(cell_metadata[,c(4,1,2,3)])
  feat_metadata <- tibble::as_tibble(SingleCellExperiment::rowData(spe),rownames = "feat_ID")
  #colnames(feat_metadata)[[1]] <- c("feat_ID")
  #rownames(raw_epxr) <- gene_metadata$gene_ID
  if (annotation) {
    rownames(raw_expr) <- c(SingleCellExperiment::rowData(spe)[, "SYMBOL"])
    #rownames(norm_expression) <- c(SingleCellExperiment::rowData(sce)[,"SYMBOL"])
  }
  gobj = Giotto::createGiottoObject(
      expression = raw_expr,
      cell_metadata = cell_metadata,
      spatial_locs = as.data.frame(SpatialExperiment::spatialCoords(spe)),
      feat_metadata = feat_metadata,
      instructions = instructions,
      dimension_reduction = GiottoClass::createDimObj(
          coordinates = SingleCellExperiment::reducedDim(spe, "reducedDim"),
          name = "PCA",
          method = "pca")
  )
  return(gobj)
}
# Convert to Giotto object
my_giotto_object <- createGiotto_fn(spe, instructions = instrs)

# Normalize
my_giotto_object <- Giotto::normalizeGiotto(my_giotto_object)

# Create nearest network
    my_giotto_object <- createNearestNetwork(
        gobject = my_giotto_object,
        type = type,
        spat_unit = "cell",
        feat_type = "rna",
        dim_reduction_to_use = "pca",
        dim_reduction_name = 'PCA',
        dimensions_to_use = 1:20, 
        k = k,
        name = 'network')

###2. Run clustering

#if(!file.exists(hmrf_folder)) dir.create(hmrf_folder, recursive = T)
message("Running ", method, " clustering")
if (method == "HMRF"){
    
    # create spatial network (required for binSpect methods)
    my_giotto_object <- Giotto::createSpatialNetwork(
        gobject = my_giotto_object,
        minimum_k = 10)
    
    # identify genes with a spatial coherent expression profile
    #km_spatialgenes <- Giotto::binSpect(my_giotto_object, bin_method = 'rank')

    #my_spatial_genes <- km_spatialgenes[1:100]$feats
    HMRF_spatial_genes <- Giotto::doHMRF(
        gobject = my_giotto_object,
        spat_unit = "cell",
        feat_type = "rna",
        betas = c(0, 2, config$beta),
        expression_values = "normalized",
        spatial_genes = rownames(spe), #my_spatial_genes,
        #dim_reduction_to_use = "pca",
        #dim_reduction_name = "PCA",
        k = n_clusters,
        name = method,
        seed = seed
        )
    #viewHMRFresults2D(gobject = my_giotto_object,
    #                    HMRFoutput = HMRF_spatial_genes,
    #                   k = 9, betas_to_view = i,
    #                  point_size = 2)
    
    
    ## Write output
    my_giotto_object <- addHMRF(
        gobject = my_giotto_object,
        HMRFoutput = HMRF_spatial_genes,
        k = k, betas_to_add = c(config$beta),
        hmrf_name = method)
    
} else if (method == "leiden"){

    my_giotto_object <- doLeidenCluster(
      my_giotto_object,
      spat_unit = "cell",
      feat_type = "rna",
      name = method,
      nn_network_to_use = type,
      network_name = "network",
      python_path = python_path,
      resolution = config$resolution,
      weight_col = "weight",
      partition_type = "RBConfigurationVertexPartition",#, "ModularityVertexPartition"),
      init_membership = NULL,
      n_iterations = 1000,
      return_gobject = TRUE,
      set_seed = TRUE,
      seed_number = seed
    )
    } else if (method == "louvain"){
    message("Version ", config$version)
    
        my_giotto_object <- doLouvainCluster(
            my_giotto_object,
            spat_unit = "cell",
            feat_type = "rna",
            version = config$version,
            name = method,
            nn_network_to_use = type,
            network_name = "network",
            python_path = python_path,
            resolution = config$resolution,
            weight_col = "weight",
            gamma = 1,
            omega = 1,
            louv_random = FALSE,
            return_gobject = TRUE,
            set_seed = TRUE,
            seed_number = seed
        )
    } else if (method == "randomwalk"){
        my_giotto_object <- doRandomWalkCluster(
            my_giotto_object,
            name = method,
            nn_network_to_use = type,
            network_name = "network",
            walk_steps = 4,
            walk_clusters = n_clusters,
            walk_weights = NA,
            return_gobject = TRUE,
            set_seed = TRUE,
            seed_number = seed
        )
    } else if (method == "sNNclust" && type == "kNN"){
        my_giotto_object <- doSNNCluster(
            my_giotto_object,
            name = method,
            nn_network_to_use = type,
            network_name = "network",
            k = config$sNNK,
            eps = config$eps,
            minPts = config$minPts,
            borderPoints = TRUE,
            return_gobject = TRUE,
            set_seed = TRUE,
            seed_number = seed
        )
    } else if (method == "kmeans"){
        my_giotto_object <- doKmeans(
            my_giotto_object,
            spat_unit = "cell",
            feat_type = "rna",
            expression_values = "normalized",#c("normalized", "scaled", "custom"),
            feats_to_use = NULL,
            dim_reduction_to_use = "pca", #c("cells", "pca", "umap", "tsne"),
            dim_reduction_name = "PCA",
            dimensions_to_use = 1:10,
            distance_method = "original",#, "pearson", "spearman", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"),
            centers = n_clusters,
            iter_max = 100,
            nstart = 1000,
            algorithm = "Hartigan-Wong",
            name = "kmeans",
            return_gobject = TRUE,
            set_seed = TRUE,
            seed_number = seed
        )
    }else if (method == "hierarchical"){
        my_giotto_object <- doHclust(
            my_giotto_object,
            spat_unit = "cell",
            feat_type = "rna",
            expression_values = "normalized", #c("normalized", "scaled", "custom"),
            feats_to_use = NULL,
            dim_reduction_to_use = "pca", #c("cells", "pca", "umap", "tsne"),
            dim_reduction_name = "PCA",
            dimensions_to_use = 1:10,
            distance_method = "pearson", #c("pearson", "spearman", "original", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"),
            agglomeration_method = "ward.D2", #c("ward.D2", "ward.D", "single", "complete", "average",    "mcquitty", "median", "centroid"),
            k = n_clusters,
            h = NULL,
            name = method,
            return_gobject = TRUE,
            set_seed = TRUE,
            seed_number = seed
        )
    }

## Write output
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
label_df <- as.data.frame(Giotto::getCellMetadata(my_giotto_object, output = "data.table"))
label_df <- data.frame(label = label_df[length(colnames(label_df))], row.names = label_df[[1]])
colnames(label_df) <- "label"

print(table(label_df$label))

write.table(label_df, file = label_file, sep = "\t", col.names = NA, quote = FALSE)


