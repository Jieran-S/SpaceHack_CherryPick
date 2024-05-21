library(tidyverse)
library(data.table)
library(RColorBrewer)
library(funkyheatmap)
library(scales)


############## Load Data -------------------------------------------------------

# Define essential path
metric_list <- c("adjusted_rand_score", "NMI", "fowlkes_mallows_score",
  "completeness_score", "homogeneity_score","v_measure_score")

continuity_list <- c("CHAOS", "PAS")

metric_name <- c("ARI", "NMI", "FMI", "COM", "HOM", "VMS") %>% set_names(metric_list)

result_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Cluster_Benchmark/Visualization/011_metric_results"
metaData    <- fread("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Cluster_Benchmark/Visualization/011_metric_results/samples_meta.tsv",header = TRUE)
data_meta   <- fread("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Cluster_Benchmark/Visualization/011_metric_results/data_meta.tsv")[-7,]
scalability <- fread("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Cluster_Benchmark/Visualization/040_Scability/scalability.tsv")[, -1]
usability   <- fread("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Cluster_Benchmark/Visualization/050_Usability/Usability_scores.csv")

method_meta <- c( 
  BayesSpace= "Bayesian", STAGATE= "GNN", BANKSY="Neighbor", spaGCN="GNN", 
  DRSC= "Bayesian", GraphST="GNN", SEDR="GNN", "scanpy-leiden"="Non-Spatial",
  "scanpy-louvain" = "Non-Spatial", "seurat-leiden"="Non-Spatial", "seurat-louvain"="Non-Spatial",
  CellCharter="Neighbor", SpaceFlow="GNN", stardust="Neighbor", bass="Bayesian", 
  precast="Bayesian")

metric_results <- lapply(metric_list, function(metric){
  # Load data from table
  data = fread(file.path(result_path, sprintf("%s.tsv", metric)), header = TRUE)
  data = data.table::melt(data, id.vars = "V1", 
                          variable.name = "sample", value.name = metric)
  setnames(data, "V1", "method_instance")
  setkey(data, sample, method_instance)
  data <- data[!is.na(get(metric)),]
  return(data)
})

contunity_results <- lapply(continuity_list, function(metric){
  # Load data from table
  data = fread(file.path(result_path, sprintf("%s.tsv", metric)), header = TRUE)
  data = data.table::melt(data, id.vars = "V1", 
                          variable.name = "sample", value.name = metric)
  setnames(data, "V1", "method_instance")
  setkey(data, sample, method_instance)
  data <- data[!is.na(get(metric)),]
  return(data)
})

############## Data manipulation -----------------------------------------------
# Get all results, combining together
data <- Reduce(function(x, y) merge(x, y, by = c("sample", "method_instance"), all = TRUE), c(metric_results, contunity_results))

# Adding metaData: Dataset, tissue, technology and methods
data[, data_github := metaData[match(data$sample, V1), dataset]
             ][, dataset := data_meta[match(data$data_github, `Name (Github)`), Name]
               ][, technology := metaData[match(data$sample, V1), technology]
                 ][, tissue := data_meta[match(data$dataset, Name), Tissue]
                   ][, method := tstrsplit(method_instance, "_", fixed=TRUE, keep=1)]

setnames(data, metric_list, metric_name)

# calculating mean across all samples
data_mean_tech <- data[, lapply(.SD, function(x) mean(x, na.rm = TRUE)), 
                  by = .(method, technology), .SDcols = metric_name] %>% 
  melt(id.vars = c("method", "technology"), variable.name = "Metrics") %>% 
  dcast(method ~ Metrics + technology)
 
data_mean_tissue <- data[, lapply(.SD, function(x) mean(x, na.rm = TRUE)), 
                         by = .(method, tissue), .SDcols = metric_name] %>% 
  melt(id.vars = c("method", "tissue"), variable.name = "Metrics") %>% 
  dcast(method ~ Metrics + tissue)

data_continuity <- data[, technology_type := ifelse(technology == "Visium", "NGS", "Imaging")
                        ][, lapply(.SD, function(x) mean(x, na.rm = TRUE)), 
                          by = .(method, technology_type), .SDcols = continuity_list] %>% 
  melt(id.vars = c("method", "technology_type"), variable.name = "Metrics")

data_continuity <- data_continuity%>% 
  dcast(method ~ Metrics + technology_type)

data_continuity_str <- data_continuity[, lapply(.SD, function(x){x = as.character(round(x, digits = 3))}), 
                                       by=method, .SDcols = colnames(data_continuity)[-1]]

setnames(data_continuity_str, colnames(data_continuity_str)[-1], paste0(colnames(data_continuity_str)[-1], "_str"))
data_continuity_str <- as.data.frame(data_continuity_str)
rownames(data_continuity_str) <- data_continuity_str$method
#data_var <- data[, lapply(.SD, function(x) var(x, na.rm = TRUE)), 
#                  by = .(method), .SDcols = metric_name]

#data_mean[, mean_var:= data_var[, rowMeans(.SD), .SDcols = metric_name]]


# Add scability scores 
scalability <- scalability[, .(runtime_Scalability = log2(mean(`runtime(min)`)+1), 
                               memory_Scalability = log2(mean(`memory(GB)`))+1), by = method]

scalability <- scalability[, runtime_Scalability := runtime_Scalability/max(runtime_Scalability)
                           ][, memory_Scalability := memory_Scalability/max(memory_Scalability)]
# Add usability scores
usability[, "scanpy-leiden" := scanpy]
usability[, "scanpy-louvain" := scanpy]
usability[, "seurat-louvain" := seurat]
usability[, "seurat-leiden" := seurat]

usability[, V1 := factor(V1, levels=unique(V1))]

usability <- usability[, lapply(.SD, mean), by = V1, .SDcols = unique(data$method)] %>% 
  melt(measure.vars =unique(data$method),  variable.name = "method", value.name = "usability") %>% 
  dcast(method ~ V1)

data_perMethod <- Reduce(function(x, y) data.table::merge.data.table(x, y, by = "method", all = TRUE), 
                         list(data_mean_tech, data_mean_tissue, data_continuity, scalability, usability))

############## Set up plotting metadata -----------------------------------------------
data_toPlot <- as.data.frame(data_perMethod)
rownames(data_toPlot) <- data_perMethod$method

if (FALSE){
  # Generate group label
  overall_group <- c( "Method",
                      rep("Accuracy:Technology", ncol(data_mean_tech)-1), 
                      rep("Accuracy:Tissue", ncol(data_mean_tissue)-1),
                      rep("Continuity", ncol(data_continuity)-1),
                      rep("Scability", ncol(scalability)-1), 
                      rep("Usability", ncol(usability)-1))
  
  geoms <- c("text", "funkyrect", "funkyrect","funkyrect", "bar", "funkyrect") %>% set_names(unique(overall_group))
  palette_list <- c("overall", "benchmark", "scaling", "stability", "overall", "qc") %>% set_names(unique(overall_group))
  
  # setnames(data_perMethod, "method", "id")
  group_level_2 <- tstrsplit(names(data_toPlot), "_", fixed = TRUE, keep = 2)[[1]][-1]
  groups <- c("Method",
              paste0("Accuracy_Technology_", group_level_2[1:(ncol(data_mean_tech)-1)]), 
              paste0("Accuracy_Tissue_", group_level_2[1:(ncol(data_mean_tissue)-1) + (ncol(data_mean_tech)-1)]),
              paste0("Continuity_", group_level_2[1: (ncol(data_continuity)-1) + (ncol(data_mean_tissue) + ncol(data_mean_tech)-2)]), 
              rep("Scability", ncol(scalability)-1), 
              rep("Usability", ncol(usability)-1))
  
  label_names <- c(" ", tstrsplit(names(data_toPlot), "_", fixed = TRUE, keep = 1)[[1]][-1])
  
  # Generating column information
  col_info <- data.table(id = names(data_toPlot), 
                         name = label_names, 
                         group = groups,
                         geom = geoms[overall_group], 
                         palette = palette_list[overall_group])
  
  group_lab <- c(" ",
                 group_level_2[1:(ncol(data_mean_tech)+ncol(data_mean_tissue)+ncol(data_continuity)-3)], 
                 rep("resources", ncol(scalability)-1), 
                 rep("Quality of code & paper",ncol(usability)-1))
  
  
  col_groups <- data.table(group = groups, 
                           level1 = overall_group, 
                           lavel2 = group_lab, 
                           palette = palette_list[overall_group])
  
  ## Generate row group and row information
  Group_Info <- c("Bayesian-based methods", 
                  "Graph Neural Network-based methods", 
                  "Neighborhood engineering methods", 
                  "Non-spatial methods") %>% set_names(unique(method_meta))
  
  row_info <- data.table(id = rownames(data_toPlot),
                         group = method_meta[rownames(data_toPlot)])
  
  row_group <- data.table(group = method_meta[rownames(data_toPlot)], 
                          Group = Group_Info[method_meta[rownames(data_toPlot)]])
  
}

############## Methods Column Information ################
colInfo_method <- data.table(id = colnames(data_toPlot)[1],
                             name = "",
                             group = "Method", 
                             geom = "text", 
                             palette = "overall", 
                             options = list(list(scale=FALSE, width=4)))

colGroup_method <- data.table(group = "Method",
                            level1 = "Method",
                            level2 = "",
                            palette = "overall")

############## Accuracy:Technology Column Information ################

col_accTec <- names(data_mean_tech)[-1]

lev2_Tec <- tstrsplit(col_accTec, "_", fixed=TRUE, keep = 2)[[1]]

colInfo_tech <- data.table(id = col_accTec,
                             name = tstrsplit(col_accTec, "_", fixed = TRUE, keep = 1)[[1]],
                             group = paste0("Accuracy_Technology_", lev2_Tec),
                             geom = "funkyrect",
                             palette = "benchmark",
                             options = list()
                             )

colGroup_tech <- data.table(group = paste0("Accuracy_Technology_", lev2_Tec),
                            level1 = "Accuracy:Technology",
                            level2 = lev2_Tec,
                            palette = "benchmark")

############## Accuracy:Tissue Column Information ################

col_accTis <- names(data_mean_tissue)[-1]

lev2_Tis <- tstrsplit(col_accTis, "_", fixed=TRUE, keep = 2)[[1]]

colInfo_tissue <- data.table(id = col_accTis,
                             name = tstrsplit(col_accTis, "_", fixed = TRUE, keep = 1)[[1]],
                             group = paste0("Accuracy_Tissue_", lev2_Tis),
                             geom = "funkyrect",
                             palette = "scaling",
                             options = list()
                             )

colGroup_tissue <- data.table(group = paste0("Accuracy_Tissue_", lev2_Tis),
                            level1 = "Accuracy:Tissue",
                            level2 = lev2_Tis,
                            palette = "scaling")

############## Continuity Column Information ################

col_cont <- rep(names(data_continuity)[-1], each = 2)
option_list <- list()
for (col in unique(col_cont)){
  option_list[[length(option_list)+1]] <- list(scale=FALSE, legend=FALSE)
  option_list[[length(option_list)+1]] <- list(label = paste0(col, "_str"), 
                                               overlay = TRUE, 
                                               size = 3, 
                                               scale = FALSE)
}

lev2_Cont <- tstrsplit(col_cont, "_", fixed=TRUE, keep = 2)[[1]]

colInfo_cont <- data.table(id = col_cont,
                             name = tstrsplit(col_cont, "_", fixed = TRUE, keep = 1)[[1]],
                             group = paste0("Continuity_",lev2_Cont),
                             geom = rep(c("rect", "text"), length(col_cont)/2),
                             palette = rep(c("stability", "white6black4"), length(col_cont)/2),
                             options = option_list # need to be defined by hand
                           )

colGroup_cont <- data.table(group = paste0("Continuity_",lev2_Cont),
                              level1 = "Continuity",
                              level2 = lev2_Cont,
                              palette = "stability")

############## Scalability Column Information ################

col_scal <- names(scalability)[-1]

lev2_scal <- rep("resources", length(col_scal))

colInfo_scal <- data.table(id = col_scal,
                           name = tstrsplit(col_scal, "_", fixed = TRUE, keep = 1)[[1]],
                           group = "Scalability",
                           geom =  "bar",
                           palette = "overall",
                           options = list(list(width=4, legend=FALSE, scale=FALSE)) # need to be defined by hand
)

colGroup_scal <- data.table(group = "Scalability",
                            level1 = "Scalability",
                            level2 = lev2_scal,
                            palette = "overall")

############## Scalability Column Information ################

col_usab <- names(usability)[-1]

colInfo_usab <- data.table(id = col_usab,
                           name = col_usab,
                           group = "Usability",
                           geom =  "funkyrect",
                           palette = "qc",
                           options = list(list(legend=FALSE, scale=FALSE)) # need to be defined by hand
)


colGroup_usab <- data.table(group = "Usability",
                            level1 = rep("Usability", length(col_usab)),
                            level2 = "Quality of code & paper",
                            palette = "qc")

############## Merge everything together ################

column_information <- rbindlist(list(colInfo_method, colInfo_tech, colInfo_tissue, colInfo_cont, colInfo_scal, colInfo_usab))
column_group <- rbindlist(list(colGroup_method, colGroup_tech, colGroup_tissue, colGroup_cont, colGroup_scal, colGroup_usab))

palettes <- dynbenchmark_data$palettes

############## Ordering data frame based on matedata ---------------------------
# Sort columns first
data_toPlot <- t(data_toPlot[, -1])
level_1_vec <- c(column_group$level1[1:(ncol(data_mean_tech)+ncol(data_mean_tissue)-2)+1], # All accuracy
                 rep(c("Continuity_bc", "Continuity_in-situ"), 2), # All continuity
                 column_group$level1[(nrow(column_group)-ncol(usability)-ncol(scalability)+3):nrow(column_group)] # ALl the rest
                 )
level_2_vec <- c(column_group$level2[1:(ncol(data_mean_tech)+ncol(data_mean_tissue)-2)+1], # All accuracy
                 rep(c("NGS", "Imaging"), 2), # All continuity
                 column_group$level2[(nrow(column_group)-ncol(usability)-ncol(scalability)+3):nrow(column_group)] # ALl the rest
)
level_3_vec <- c(column_information$name[1:(ncol(data_mean_tech)+ncol(data_mean_tissue)-2)+1], # All accuracy
                 rep(c("CHAOS", "PAS"), each=2), # All continuity
                 column_information$name[(nrow(column_group)-ncol(usability)-ncol(scalability)+3):nrow(column_group)] # ALl the rest
)

level_1 <- factor(level_1_vec, levels = unique(level_1_vec))
level_2 <- factor(level_2_vec, levels = c("Visium", "STARmap+", "MERFISH", "Xenium",
                                          unique(level_2_vec)[5:length(unique(level_2_vec))]))
level_3 <- factor(level_3_vec, levels = unique(level_3_vec))

data_toPlot <- data_toPlot[order(level_1, level_2, level_3), ]

# Then the rows
data_toPlot <- t(data_toPlot)
# row_1 <- factor(row_group$Group, levels = Group_Info)
# Rank rows based on average ARI
metMeans <- data[, lapply(.SD, mean), by = method, .SDcols = metric_name]
metRanks <- metMeans[, lapply(.SD, frank), .SDcols = metric_name]
rowRanks <- - rowMeans(metRanks) %>% set_names(metMeans$method)

# rowMean <- -rowMeans(data_toPlot[, 1:(ncol(data_mean_tech)+ncol(data_mean_tissue)-2)], na.rm = TRUE)
data_toPlot <- data_toPlot[order(rowRanks),]

data_toPlot <- as.data.frame(data_toPlot)
data_toPlot <- cbind(method = rownames(data_toPlot), data_toPlot)

match_colnames <- lapply(match(colnames(data_toPlot), column_information$id), 
                         function(x){
                           if (x > ncol(data_mean_tech)+ncol(data_mean_tissue)-1 &
                               x < ncol(data_mean_tech)+ncol(data_mean_tissue)+8){
                             return(c(x, x+1))
                           } else {
                             return(x)
                           }
                         }) %>% unlist()

data_toPlot_merge <- data_toPlot %>% right_join(data_continuity_str)
rownames(data_toPlot_merge) <- rownames(data_toPlot)

g <- funky_heatmap(
  data = data_toPlot_merge,
  column_info = column_information[match_colnames,],
  column_groups = column_group[match_colnames,],
  #row_info = row_info[match(rownames(data_toPlot), row_info$id),],
  #row_groups = row_group[match(rownames(data_toPlot), row_info$id),],
  palettes = palettes,
  scale_column = TRUE,
  position_args = position_arguments(col_annot_offset = 2.5, col_annot_angle = 30)
)

g
