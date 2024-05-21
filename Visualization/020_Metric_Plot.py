#%% Define packages
import os
from pathlib import Path
import pandas as pd
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from PyComplexHeatmap import *

#%% Define get ARI function
def get_method_metric(methods, data_path, samples, metric = "ARI"):
    import os
    import json
    import pandas as pd
    
    def extract_result(dir, method_instance):
        with open(dir, "r") as f:
            result = json.load(f)
        if isinstance(result, list):
            result = pd.DataFrame(result, columns=["domain", method_instance])
            result.set_index("domain", inplace=True)
            result = result.iloc[:,0]
        return result

    results_dict = {}
    for sample in samples:
        sample_result = {}
        avail_methods = set(methods).intersection(os.listdir(os.path.join(data_path, sample)))
        if len(avail_methods) == 0:
            continue
        else:
            for method in avail_methods:
                method_dir = os.path.join(data_path, sample, method)
                if metric in os.listdir(method_dir):
                    sample_result[method] = extract_result(os.path.join(method_dir, metric, "results.txt"), 
                                                 method)
                else:
                    configs = [x for x in os.listdir(method_dir) if x.startswith("config")]
                    for config in configs:
                        if metric in os.listdir(os.path.join(method_dir, config)):
                            config_name = config.split("config_")[1]
                            if "leiden" in config:
                                mth_ins = f"{method}-leiden_{config_name}"
                            elif "louvain" in config:
                                mth_ins = f"{method}-louvain_{config_name}"
                            else:
                                mth_ins = f"{method}_{config_name}"
                            sample_result[mth_ins] = extract_result(os.path.join(method_dir, config, metric, "results.txt"), 
                                                          mth_ins)
        results_dict[sample] = pd.DataFrame.from_dict(sample_result, orient="index")
    ari_df = pd.concat(results_dict, axis=1)
    ari_df = ari_df.sort_index(axis=0)
    return ari_df

#%% Define methods and datasets
allMethods =["BayesSpace", 
              "STAGATE", 
              "BANKSY", 
              "spaGCN", 
              "DRSC", 
              "GraphST", 
              "SEDR", 
              "scanpy", 
              "seurat", 
              "CellCharter",
              "SpaceFlow",
              "stardust",
              "bass", 
              "precast"]

allMetric = [
    "ARI",
    "domain-specific-f1",
    "CHAOS",                            #phyiscal coord only
    "PAS",      # buggy 
    "Completeness",
    "Entropy",
    "FMI",
    "Homogeneity",
    "MCC",
    "NMI",
    "jaccard",
    "V_measure",                        #Config, GT
    "LISI",                             #Config, embed, GT
    "cluster-specific-silhouette",      #embed, no GT
    "Calinski-Harabasz",
    "Davies-Bouldin",
]

metric = "domain-specific-f1"

# Define directories
dataset_path = Path("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Cluster_Benchmark/data/")
save_dir = Path(f"/work/PRTNR/CHUV/DIR/rgottar1/spatial/Cluster_Benchmark/Visualization/020_{metric}")
save_dir.mkdir(exist_ok=True)


# data_name = "libd_dlpfc"
data_name = "visium_breast_cancer_SEDR"
# data_name = "visium_chicken_heart"
# data_name = "STARmap_plus"  
# data_name = "abc_atlas_wmb_thalamus"
# data_name = "xenium-mouse-brain-SergioSalas"
# data_name = "stereoseq_developing_Drosophila_embryos_larvae"

data_list = [
    "libd_dlpfc", "visium_breast_cancer_SEDR", "visium_chicken_heart", 
    "STARmap_plus", "abc_atlas_wmb_thalamus", "xenium-mouse-brain-SergioSalas",
]


data_path = dataset_path / data_name

wes_palette = [
    "#D8A499",  "#5FBFF9", "#A4243B", "#E7BB41",
    "#2F343B", "#FAD77B", "#7D4E57", "#3E442B",
    "#87B2C7", "#9B5094", "#6A994E", "#F1BB7B",
    "#FD6467", "#DF9D5A", "#F2D3BC", "#DAD6D6",
]

#%%
result_list = []
for data in data_list: 
    data_path = dataset_path / data
    samples = [m for m in os.listdir(data_path) if os.path.isdir(os.path.join(data_path, m)) and not m.startswith(".")]
    results = get_method_metric(methods= allMethods,
                            data_path=data_path,
                            samples=samples,
                            metric=metric)
    result_list.append(results)

result = result_list[0]
for result_df in result_list[1:]:
    result = result.join(result_df, how="outer")

result = result.droplevel(1, axis = 1)
result.to_csv(f"/work/PRTNR/CHUV/DIR/rgottar1/spatial/Cluster_Benchmark/Visualization/011_metric_results/{metric}.tsv", sep="\t", index_label="")

#%% Generate results
if os.path.exists(save_dir / f"{metric}_{data_name}.tsv"):
    results =  pd.read_table(save_dir / f"{metric}_{data_name}.tsv")
else:
    samples = [m for m in os.listdir(data_path) if os.path.isdir(os.path.join(data_path, m)) and not m.startswith(".")]
    results = get_method_metric(methods= allMethods,
                            data_path=data_path,
                            samples=samples,
                            metric=metric)

    # if results.shape[1] > 1:
    results = results.dropna(axis=1, how='all').T
    index_value = results.index.get_level_values(1) if isinstance(results.index, pd.MultiIndex) else results.index
    results = results[index_value.notna()]
    #results = results.T.drop_duplicates(keep="last").T

    results.to_csv(save_dir / f"{metric}_{data_name}.tsv", sep="\t", index_label=False)

# results.fillna(0.0, inplace=True)
index_value = results.index.get_level_values(1) if isinstance(results.index, pd.MultiIndex) else results.index
results = results[index_value.notna()].astype(float)
results = results[~index_value.str.contains('unassigned')].astype(float)
#%% Boxplots
id_var = ["samples", "domains"] if isinstance(results.index, pd.MultiIndex) else "samples"
datall = pd.melt(results.reset_index(names=id_var), id_vars=id_var, var_name='method_instance', value_name=metric)
datall["method"] = datall['method_instance'].str.split('/').str[0]
datall = datall.dropna(axis=0)
datall["plot_sample"] = datall["samples"] if isinstance(results.index, pd.MultiIndex) else "all"

for sample_id in datall["plot_sample"].unique():
    data = datall[datall["samples"] == sample_id]
    plt.figure(figsize=(20, 9))  # Adjust figure size if needed
    ax = sns.boxplot(data=data,  x = "method_instance", y = metric, hue = "method", dodge=False, palette = wes_palette)
    sns.stripplot(data=data, x = "method_instance", y = metric, ax=ax, hue = "method",
                jitter=True, size=5, dodge=False, palette = wes_palette)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:], labels[:], title='method', bbox_to_anchor=(1, 1.02), loc='upper left')
    plt.title(f"{data_name}: {metric}")
    plt.xticks(rotation=90)  # Rotate x-axis labels for better readability
    # plt.show()
    plt.savefig(save_dir / f"{metric}_{data_name}_{sample_id}.png", bbox_inches="tight")

#%% Barplot
id_var = ["samples", "domains"] if isinstance(results.index, pd.MultiIndex) else "samples"
datall = pd.melt(results.reset_index(names=id_var), id_vars=id_var, var_name='method_instance', value_name=metric)
datall["method"] = datall['method_instance'].str.split('/').str[0]
datall = datall.dropna(axis=0)
datall["plot_sample"] = datall["samples"] if isinstance(results.index, pd.MultiIndex) else "all"

for sample_id in datall["plot_sample"].unique():
    data = datall[datall["samples"] == sample_id]
    plt.figure(figsize=(25, 9))  # Adjust figure size if needed
    ax = sns.barplot(data=data,  x = "method_instance", y = metric, hue = "method", dodge=False, palette = wes_palette)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:], labels[:], title='method', bbox_to_anchor=(1, 1.02), loc='upper left')

    # ax.legend(handles[:], labels[:], title='method', bbox_to_anchor=(1, 1.02), loc='upper left')
    plt.title(f"{data_name}: {metric}")
    plt.xticks(rotation=90)  # Rotate x-axis labels for better readability
    # plt.show()
    plt.savefig(save_dir / f"{metric}_Bar_{data_name}_{sample_id}.png", bbox_inches="tight")



#%% Heatmap
col_anno = pd.DataFrame({"method": results.columns.str.split('/').str[0], 
                         # "mean": results.label
                         }, 
                         index=results.columns)    # col_dic = dict(zip(colored_col.unique(), wes_palette[:len(colored_col.unique())]))


if isinstance(results.index, pd.MultiIndex):
    domain_size = []
    for index in results.index:
        labels = pd.read_table(data_path / index[0] / "labels.tsv")["label"]
        domain_size.append(np.sum(labels == index[1]))

    row_anno = pd.DataFrame({
        "samples" : results.index.get_level_values(0),
        "domains" : results.index.get_level_values(1),
        "size" : domain_size,
    }, index= results.index)

    # Group the methods by SAMPLE (level0) or DOMAIN (level1), row_split changes correspondingly
    col_boxplot = results.groupby(level=1).mean().T
    row_split = row_anno.domains
    
    box = anno_boxplot(results, cmap="cividis", height=20, legend=False)
    box.ylim = box.ylim[::-1]
    row_ha = HeatmapAnnotation(
        size= anno_barplot(row_anno[["size"]], height=10, legend=True),
        domain_res=box, #anno_boxplot(results, cmap="cividis", height=20, legend=False),
        samples=anno_simple(row_anno.samples),
        domains=anno_simple(row_anno.domains, legend=False), 
        sample=anno_label(row_split, merge=True, rotation=0), 
        axis=0, orientation='right', label_side="bottom", hgap=0.3)

else:
    row_anno = None
    row_split = None
    col_boxplot = results.T
    row_ha = HeatmapAnnotation(all_method_result=anno_boxplot(results, cmap="cividis", height=30, legend=False),
                               axis=0)

plt.figure(figsize=(15, 10))
col_ha = HeatmapAnnotation(method=anno_label(col_anno.method, merge=True,rotation=30),
                            methods=anno_simple(col_anno.method, legend=False),
                            Domain_metric=anno_boxplot(col_boxplot, cmap = "plasma", height=30, legend=False),
                            axis=1,)
cm = ClusterMapPlotter(data=results, 
top_annotation=col_ha, 
right_annotation=row_ha,
col_split = col_anno.method, 
col_cluster =False, 
row_split=row_split, 
row_cluster=False, 
col_split_gap = 1, 
show_rownames=False, 
# row_names_side="left",
label = metric,
# legend_side="left",
# xticklabels = True,
cmap = "viridis")

plt.suptitle(f"{data_name}: {metric}", y = 0.98)
plt.savefig(save_dir / f"{metric}_Hmp_{data_name}_bydomain.png", bbox_inches="tight")

# %% Change the plot
results = results.dropna(axis=1, how = "any")

results_trans = results.transpose()

col_anno = pd.DataFrame({"method": results.columns.str.split('_').str[0], 
                         # "mean": results.label
                         }, 
                         index=results.columns)    # col_dic = dict(zip(colored_col.unique(), wes_palette[:len(colored_col.unique())]))


if isinstance(results.index, pd.MultiIndex):
    domain_size = []
    for index in results.index:
        labels = pd.read_table(data_path / index[0] / "labels.tsv")["label"]
        domain_size.append(np.sum(labels == index[1]))

    row_anno = pd.DataFrame({
        "samples" : results.index.get_level_values(0),
        "domains" : results.index.get_level_values(1),
        "size" : domain_size,
    }, index= results.index)

    # Group the methods by SAMPLE (level0) or DOMAIN (level1), row_split changes correspondingly
    col_boxplot = results.groupby(level=1).mean().T
    row_split = row_anno.domains
    
    box = anno_boxplot(results, cmap="winter", height=40, legend=False)
    box.ylim = box.ylim[::-1]
    col_ha = HeatmapAnnotation(
        sample=anno_label(row_split, merge=True, rotation=30), 
        domains=anno_simple(row_anno.domains, legend=False), 
        no_cells= anno_barplot(row_anno[["size"]], height=40, legend=False, cmap="Wistia"),
        # domain_F1=box, #anno_boxplot(results, cmap="cividis", height=20, legend=False),
        # samples=anno_simple(row_anno.samples),
        axis=1, vgap=0.3)

else:
    row_anno = None
    row_split = None
    col_boxplot = results.T
    row_ha = HeatmapAnnotation(all_method_result=anno_boxplot(results, cmap="cividis", height=30, legend=False),
                               axis=0)

# Define column order
col_order = results.groupby(level=1).mean().mean(axis=1).sort_values().index.tolist()

plt.figure(figsize=(9, 6))
row_ha = HeatmapAnnotation(
    method=anno_label(col_anno.method, merge=True,rotation=0, fontsize = 7),
    m=anno_simple(col_anno.method, legend=False),
    # Domain_metric=anno_boxplot(col_boxplot, cmap = "plasma", height=30, legend=False),
    axis=0, label_side="bottom")
cm = ClusterMapPlotter(data=results_trans, 
top_annotation=col_ha, 
left_annotation=row_ha,
row_split = col_anno.method, 
col_cluster =False, 
col_split=row_split, 
row_cluster=False, 
col_split_order=col_order,
#col_split_gap = 1, 
row_split_gap=0.15,
show_rownames=False, 
label = "Domain F1",
# legend_side="left",
# xticklabels = True,
cmap = "viridis")

#plt.suptitle(f"{data_name}: {metric}", y = 0.98)
plt.savefig(save_dir / f"{metric}_Hmp_{data_name}_bydomain.png", bbox_inches="tight")

# %% Spearman correlation
if isinstance(results.index, pd.MultiIndex):
    domain_size = []
    for index in results.index:
        labels = pd.read_table(data_path / index[0] / "labels.tsv")["label"]
        domain_size.append(np.sum(labels == index[1]))

results_nonSpatial  = pd.DataFrame()
columns = ["scanpy-leiden", "scanpy-louvain", "seurat-leiden", "seurat-louvain"]
for col in columns:
    columns_to_keep = results.filter(regex=col).columns
    print(columns_to_keep)
    results_nonSpatial[col] = results[columns_to_keep].mean(axis=1)

results_nonSpatial["n_cells"] = domain_size

corr_df = pd.melt(results_nonSpatial, id_vars="n_cells")

from scipy.stats import spearmanr

coef, p = spearmanr(corr_df["value"], corr_df["n_cells"])

sns.set(style="whitegrid")
plt.figure(figsize=(5.5, 5))
sns.scatterplot(x='value', y='n_cells', hue='variable', data=corr_df, s=100)


# Annotate the plot with the Spearman correlation coefficient and p-value
#plt.annotate(f'Spearman correlation: {coef:.2f}\nP-value: {p:.2e}', 
#             xy=(0.05, 0.95), xycoords='axes fraction', fontsize=12,
#             horizontalalignment='left', verticalalignment='top')

# Set plot labels and title
plt.xlabel('')
plt.ylabel('Number of spots in domain')
# handles, labels = plt.get_legend_handles_labels()
plt.legend(title='method', bbox_to_anchor=(1, 1.02), loc='upper left')

# plt.title('Correlation between Mean F1 for non-spatial methods and Number of Cells')

plt.show()

# %% Multiple-test and p value correction
from statsmodels.stats.multitest import multipletests

df = results_nonSpatial

# Calculate Spearman correlations and p-values for each metric with 'y'
resultspear = {'method': [], 'pearson correlation': [], 'p value': []}
for column in df.columns[:-1]:  
    corr, p_value = spearmanr(df[column], df['n_cells'])
    resultspear['method'].append(column)
    resultspear['pearson correlation'].append(corr)
    resultspear['p value'].append(p_value)

results_df = pd.DataFrame(resultspear)

# Apply Bonferroni correction
alpha = 0.05
results_df['bonferroni corrected p value'] = multipletests(results_df['p value'], alpha=alpha, method='bonferroni')[1]

# Apply FDR correction (Benjamini-Hochberg)
results_df['fdr corrected p value'] = multipletests(results_df['p value'], alpha=alpha, method='fdr_bh')[1]

fig, ax = plt.subplots(figsize=(6, 2))  # Set the figure size
ax.axis('tight')
ax.axis('off')
table = ax.table(cellText=results_df.values, colLabels=results_df.columns, cellLoc='center', loc='center')

plt.show()

# %%
