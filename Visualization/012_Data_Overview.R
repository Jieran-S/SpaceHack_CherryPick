library(tidyverse)
library(data.table)
library(RColorBrewer)
library(funkyheatmap)
library(kableExtra)

######### Read Data in ---------------------------------------------------------

data <- fread("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Cluster_Benchmark/Visualization/011_metric_results/data_meta.tsv")[-7,]

wes_palette = c(
  "#A4243B", "#F98400", "#FFDE05", "#0B775E", "#09ACEC",  
  "#2F343B", "#7D4E57", "#E7BB41",
  "#87B2C7", "#9B5094", "#F1BB7B", "#3E442B",
  "#FD6467", "#DF9D5A", "#F2D3BC", "#DAD6D6"
)


######### Adding gene, cell and nclust info ------------------------------------

sample_meta <- fread("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Cluster_Benchmark/Visualization/011_metric_results/samples_meta.tsv")
sample_mean <- sample_meta[, lapply(.SD, mean), by = dataset, .SDcols = c("n_clusters", "n_genes", "n_cells")]


data[, c("Genes", "Genes_num", 
         "Cells", "Cells_num", 
         "N_Cluster", "N_Cluster_num") := lapply(rep(c("n_genes", "n_cells", "n_clusters"), each=2), 
                                                  function(x){sample_mean[match(data$`Name (Github)`, dataset), get(x)]})]
data = data[, -1]


data[, N_Cluster_num := as.integer(N_Cluster_num)][, Genes_num := as.integer(Genes)][, Cells_num := as.integer(Cells)]

setcolorder(data, c(colnames(data)[1:8], "Genes", "Genes_num", 
                 "Cells", "Cells_num", 
                 "N_Cluster", "N_Cluster_num"))

######## Plot funky heatmap ----------------------------------------------------

col_info <- data.table(id = colnames(data),
                       geom = c(rep("text", 8), rep(c("bar","text"), 3)),
                       name = rep("", 14),
                       group = c(colnames(data)[1:8], rep(c("Genes", "Cells", "N_Cluster"), each = 2)),
                       palette = c(rep("black", 7), "Histology", "greens", "black", "blues","black", "oranges", "black")
                       #width = c(rep(3, 8), rep(5, 3))
                       )

column_groups <- data.table(group = c(colnames(data)[1:8], rep(c("Genes", "Cells", "N_Cluster"), each = 2)), 
                            palette = "title")

palettes <- list(
  # features = c(FULL = "#4c4c4c", HVG = "#006300"),
  blues = "Blues",
  greens = "Greens",
  title = wes_palette,
  Histology = c(No = "darkred", Yes = "#006300"),
  oranges = rev(RColorBrewer::brewer.pal(9, "Oranges")),
  #greys = "Greys",
  black = rep("black", 6)
)

g <- funky_heatmap(
  data = data,
  column_info = col_info,
  palettes = palettes,
  position_args = position_arguments(
    col_annot_offset = 1
  ),
  column_groups = column_groups,
  scale_column = TRUE
)

g


ggplot(sample_meta, aes(x = dataset, y = n_cells)) + 
  geom_boxplot() +
  theme_bw()
 
