library(tidyverse)
library(data.table)
library(ggbeeswarm)
library(RColorBrewer)
library(wesanderson)
library(ggridges)


theme_pub <- function(
    font_size = 14,
    font_family = "",
    line_size = .5,
    rel_small = 12/14,
    rel_tiny = 11/14,
    rel_large = 16/14
) {
  half_line <- font_size / 2
  small_size <- rel_small * font_size
  
  # work off of theme_grey just in case some new theme element comes along
  theme_grey(base_size = font_size, base_family = font_family) %+replace%
    theme(
      line              = element_line(colour = "black", size = line_size, linetype = 1, lineend = "butt"),
      rect              = element_rect(fill = NA, colour = NA, size = line_size, linetype = 1),
      text              = element_text(family = font_family, face = "plain", colour = "black",
                                       size = font_size, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = .9,
                                       margin = margin(), debug = FALSE),
      
      axis.line         = element_line(colour = "black", size = line_size, lineend = "square"),
      axis.line.x       = NULL,
      axis.line.y       = NULL,
      axis.text         = element_text(colour = "black", size = small_size),
      axis.text.x       = element_text(margin = margin(t = small_size / 4), vjust = 1),
      axis.text.x.top   = element_text(margin = margin(b = small_size / 4), vjust = 0),
      axis.text.y       = element_text(margin = margin(r = small_size / 4), hjust = 1),
      axis.text.y.right = element_text(margin = margin(l = small_size / 4), hjust = 0),
      axis.ticks        = element_line(colour = "black", size = line_size),
      axis.ticks.length = unit(half_line / 2, "pt"),
      axis.title.x      = element_text(
        margin = margin(t = half_line / 2),
        vjust = 1
      ),
      axis.title.x.top  = element_text(
        margin = margin(b = half_line / 2),
        vjust = 0
      ),
      axis.title.y      = element_text(
        angle = 90,
        margin = margin(r = half_line / 2),
        vjust = 1
      ),
      axis.title.y.right = element_text(
        angle = -90,
        margin = margin(l = half_line / 2),
        vjust = 0
      ),
      
      
      legend.background = element_blank(),
      legend.spacing    = unit(font_size, "pt"),
      legend.spacing.x  = NULL,
      legend.spacing.y  = NULL,
      legend.margin     = margin(0, 0, 0, 0),
      legend.key        = element_blank(),
      legend.key.size   = unit(1.1 * font_size, "pt"),
      legend.key.height = NULL,
      legend.key.width  = NULL,
      legend.text       = element_text(size = rel(rel_small)),
      legend.text.align = NULL,
      legend.title      = element_text(hjust = 0),
      legend.title.align = NULL,
      legend.position   = "right",
      legend.direction  = NULL,
      legend.justification = c("left", "center"),
      legend.box        = NULL,
      legend.box.margin =  margin(0, 0, 0, 0),
      legend.box.background = element_blank(),
      legend.box.spacing = unit(font_size, "pt"),
      
      panel.background  = element_blank(),
      panel.border      = element_blank(),
      panel.grid        = element_blank(),
      panel.grid.major  = NULL,
      panel.grid.minor  = NULL,
      panel.grid.major.x = NULL,
      panel.grid.major.y = NULL,
      panel.grid.minor.x = NULL,
      panel.grid.minor.y = NULL,
      panel.spacing     = unit(half_line, "pt"),
      panel.spacing.x   = NULL,
      panel.spacing.y   = NULL,
      panel.ontop       = FALSE,
      
      strip.background  = element_rect(fill = "grey80"),
      strip.text        = element_text(
        size = rel(rel_small),
        margin = margin(half_line / 2, half_line / 2,
                        half_line / 2, half_line / 2)
      ),
      strip.text.x      = NULL,
      strip.text.y      = element_text(angle = -90),
      strip.placement   = "inside",
      strip.placement.x =  NULL,
      strip.placement.y =  NULL,
      strip.switch.pad.grid = unit(half_line / 2, "pt"),
      strip.switch.pad.wrap = unit(half_line / 2, "pt"),
      
      plot.background   = element_blank(),
      plot.title        = element_text(
        face = "bold",
        size = rel(rel_large),
        hjust = 0.5, vjust = 1,
        margin = margin(b = half_line)
      ),
      plot.subtitle     = element_text(
        size = rel(rel_small),
        hjust = 0, vjust = 1,
        margin = margin(b = half_line)
      ),
      plot.caption      = element_text(
        size = rel(rel_tiny),
        hjust = 1, vjust = 1,
        margin = margin(t = half_line)
      ),
      plot.tag           = element_text(
        face = "bold",
        hjust = 0, vjust = 0.7
      ),
      plot.tag.position = c(0, 1),
      plot.margin       = margin(half_line, half_line, half_line, half_line),
      
      complete = TRUE
    )
}


############## Load Data -------------------------------------------------------

# Define essential path
metric_list <- c(
  "adjusted_rand_score", 
  "NMI",
  "fowlkes_mallows_score",
  "completeness_score",
  "homogeneity_score",
  "v_measure_score")


result_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Cluster_Benchmark/Visualization/011_metric_results"

metric <- metric_list[2]

# Load data from table
resultData = fread(file.path(result_path, sprintf("%s.tsv", metric)), header = TRUE)
metaData = fread("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Cluster_Benchmark/Visualization/011_metric_results/samples_meta.tsv",
                       header = TRUE)
data_meta   <- fread("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Cluster_Benchmark/Visualization/011_metric_results/data_meta.tsv")[-7,]

metaData[, tissue := data_meta[match(metaData$dataset, `Name (Github)`), Tissue]]
metaData[, dataset_name := data_meta[match(metaData$dataset, `Name (Github)`), Name]]

method_meta <- c(
  BayesSpace= "Bayesian", 
  STAGATE= "GNN", 
  BANKSY="Neighbor", 
  spaGCN="GNN", 
  DRSC= "Bayesian",
  GraphST="GNN", 
  SEDR="GNN", 
  "scanpy-louvain"="Non-Spatial", 
  "seurat-louvain"="Non-Spatial",
  "scanpy-leiden"="Non-Spatial",
  "seurat-leiden"="Non-Spatial",
  CellCharter="Neighbor",
  SpaceFlow="GNN",
  stardust="Neighbor",
  bass="Bayesian", 
  precast="Bayesian"
)

############## Data manipulation -----------------------------------------------

data <- data.table::melt(resultData, id.vars = "V1", 
                         variable.name = "sample", value.name = "metric_value")
setnames(data, "V1", "method_instance")
data[, dataset := metaData[match(data$sample, V1), dataset_name]]
data[, dataset_github := metaData[match(data$sample, V1), dataset]]
data[, tissue := metaData[match(data$sample, V1), tissue]]
data[, tissue_type := ifelse(tissue == "Brain", "Brain", "Non-Brain")]
data[, technology := metaData[match(data$sample, V1), technology]]
data[, method := tstrsplit(method_instance, "_", fixed=TRUE, keep=1)]
data[, method_type := method_meta[method]]
data <- data[!is.na(metric_value), ]

data_mean <- data[, .(mean = mean(metric_value)), by = method
                  ][order(-mean), ]

method_order <- c(1:nrow(data_mean)) %>% set_names(data_mean$method)

data_region <- data[, .(max = max(metric_value), min = min(metric_value)), by = .(method, dataset)
                    ][, x_pos := method_order[method]][]

wes_palette = c(
  "#A4243B", "#F98400", "#FFDE05", "#0B775E", "#09ACEC",  
  "#2F343B", "#7D4E57", "#E7BB41",
  "#87B2C7", "#9B5094", "#F1BB7B", "#3E442B",
  "#FD6467", "#DF9D5A", "#F2D3BC", "#DAD6D6"
)

############## violin plot: dot with different datasets ------------------------

wes_dic <- wes_palette[1:length(unique(data$dataset))] %>% set_names(unique(data$dataset))

#700, 520

plot_vio_point <- ggplot(data, aes(x = factor(method, levels=data_mean$method))) +
  stat_summary(aes(y = metric_value), fun = "mean", geom = "crossbar", width = 0.5) + 
  geom_ribbon(data = data_region[dataset=="XMB"], aes(x=x_pos, ymin = min, ymax=max, fill = dataset), alpha = 0.5) +
  ggbeeswarm::geom_quasirandom(aes(y = metric_value, color = dataset), shape = 16, size = 1.5) +
  stat_summary(aes(y = metric_value), fun = "mean", geom = "crossbar", width = 0.5) + 
  scale_fill_manual(data$dataset, values = wes_dic) +
  scale_color_manual(data$dataset, values = wes_dic) +
  labs(x = "Method", y = metric, title = metric) + 
  # scale_x_discrete("", labels = method_instance) +
  # scale_y_continuous(metric, expand = c(0, 0), limits = c(0, 1)) +
  theme_pub(font_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.key.size = unit(20, "pt"),
        legend.position = "right", 
        legend.justification = "left",
        panel.grid.major.y = element_line(size = 0.5, color = "#DDDDDD")) +
  guides(color = guide_legend(ncol = 1, title="Datasets"), 
         fill = "none")

plot_vio_point


##### Another filled violin plot -----------------------------------------------

dataset_mean <- data[,.(mean = mean(metric_value)), by = dataset
                     ][order(-mean),
                       ]

bw <- 0.05
densities <-
  data %>%
  group_by(method, dataset) %>%
  summarise(density = list(density(metric_value, bw = bw, from = 0, to = 1, n = 100))) %>%
  mutate(x = map(density, "x"), y = map(density, "y")) %>%
  unnest(x, y) %>%
  ungroup()

densities_stacked <-
  densities %>%
  group_by(method, x) %>%
  mutate(dataset = factor(dataset, dataset_mean$dataset)) %>% # set order of trajectory types
  arrange(dataset) %>%
  mutate(norm = sum(y), y = y * y, y = y / sum(y) * norm, y = ifelse(is.na(y), 0, y)) %>% # normalise between 0 and 1
  mutate(ymax = cumsum(y), ymin = lag(ymax, default = 0)) %>%
  ungroup() %>%
  group_by(method) %>%
  mutate(ymin = ymin / max(ymax), ymax = ymax / max(ymax)) %>% # normalise so that the maximal density is 1
  ungroup()

densities_violin <-
  densities_stacked %>%
  group_by(method, x) %>%
  mutate(ymax_violin = ymax - max(ymax)/2, ymin_violin = ymin - max(ymax)/2) %>%
  ungroup()

densities_violin$method_id <- method_order[densities_violin$method]

plot_vio_density <-
  ggplot(densities_violin) +
  geom_ribbon(
    aes(
      x,
      ymin = ymin_violin + as.numeric(method_id),
      ymax = ymax_violin + as.numeric(method_id),
      fill = dataset,
      group = paste0(method, dataset),
      # alpha = dataset %in% !!trajectory_types_oi
    ), position = "identity"
  ) +
  scale_fill_manual(values=wes_dic) +
  geom_point(aes(y = method_order, x = mean), data = data_mean, size = 20, shape = 45, color = "black") +
  scale_y_continuous(NULL, breaks = seq_along(method_order), labels = data_mean$method, expand = c(0, 0)) +
  scale_x_continuous("Overall score", limits = c(0, 1), expand = c(0, 0)) +
  coord_flip() +
  labs(x = "Method", y = metric, title = metric) + 
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "right", legend.justification = "left") +
  guides(fill = guide_legend(ncol = 1), alpha = FALSE)

plot_vio_density

##### Ridge plot recreation -----------------------------------------------

bandwidth <- 1

# Ridge plot of the method performance separated by different method type (row) 
# and different technology (col)

data[, technology_type := ifelse(technology == "Visium", "Barcode", "in-situ")]

data_vline_mean = data[, mean(metric_value), by = technology_type]

method_order_ridge = data[, .(mean = mean(metric_value)), by = .(method, method_type)
                    ][order(method_type, -mean), method]

plot_method_acccuracy <- ggplot(data, aes(metric_value, y = factor(method, levels = method_order_ridge))) +
  ggridges::geom_density_ridges2(aes(fill = method_type)) +
  theme_pub() +
  xlab(metric) + 
  ylab("method type") + 
  scale_fill_manual(values = wes_palette[-3]) + 
  geom_vline(data = data_vline_mean, aes(xintercept = V1), color = "black", linetype = "dashed") +
  facet_wrap(~technology_type) + 
  theme(legend.position = "top", legend.justification = "center")


plot_method_acccuracy



save_path = "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Cluster_Benchmark/Visualization/022_metric_violin"
dir.create(save_path, showWarnings = FALSE)

ggsave(file.path(save_path, sprintf("%s_pnt_dataset.png", metric)), plot = plot_vio_point, height = 520/3, width = 800/3, dpi = 600, units = "mm")
ggsave(file.path(save_path, sprintf("%s_dst_dataset.png", metric)), plot = plot_vio_density, height = 520/3, width = 800/3, dpi = 600, units = "mm")

