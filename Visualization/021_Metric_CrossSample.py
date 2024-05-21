#%% Define packages
import os
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from PyComplexHeatmap import *

#%% Define get ARI function
def get_data_metric(methods, data_path, samples, metric = "ARI"):
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
                            mth_ins = method + "/" + config
                            sample_result[mth_ins] = extract_result(os.path.join(method_dir, config, metric, "results.txt"), 
                                                          mth_ins)
        results_dict[sample] = pd.DataFrame.from_dict(sample_result, orient="index")
    ari_df = pd.concat(results_dict, axis=1)
    ari_df = ari_df.sort_index(axis=0)
    return ari_df

#%% Define methods and datasets
allMethods = ["BayesSpace", 
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
              "precast",
              "stardust",
              "bass"]

allMetric = [
    "ARI",
    "domain-specific-f1",
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
    "CHAOS",                            #phyiscal coord only
    "PAS"      # buggy 
]

#metric = allMetric[0]

#metric = "adjusted_rand_score"
#metric = "completeness_score"
#metric = "fowlkes_mallows_score"
#metric = "homogeneity_score"
#metric = "normalized_mutual_info_score"
metric = "v_measure_score"

dataset_path = Path("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Cluster_Benchmark/data/")
metric_res_path = Path("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Cluster_Benchmark/Visualization/011_metric_results")
save_dir = Path(f"/work/PRTNR/CHUV/DIR/rgottar1/spatial/Cluster_Benchmark/Visualization/021_CrossSample")

data_list = [
    "libd_dlpfc",
    "visium_breast_cancer_SEDR",
    "visium_chicken_heart",
    "STARmap_plus",
    "abc_atlas_wmb_thalamus",
    "xenium-mouse-brain-SergioSalas",
    # "stereoseq_developing_Drosophila_embryos_larvae"
]
#%% Generate results
def get_method_metric(data_list, methods, metric, prior_dir=None):
    from pathlib import Path
    import os

    if prior_dir is None:
        prior_dir = Path(f"/work/PRTNR/CHUV/DIR/rgottar1/spatial/Cluster_Benchmark/Visualization/020_{metric}")

    results_dic = {}
    for data in data_list:
#        if os.path.isfile(prior_dir / f"{metric}_{data}.tsv"):
#            results = pd.read_table(prior_dir / f"{metric}_{data}.tsv", index_col=0)
#            results_list.append(results)
#        else: 
        data_path = dataset_path / data
        samples = [m for m in os.listdir(data_path) if os.path.isdir(os.path.join(data_path, m)) and not m.startswith(".")]
        results = get_data_metric(methods= methods,
                        data_path=data_path,
                        samples=samples,
                        metric=metric)

        # if results.shape[1] > 1:
        results = results.dropna(axis=1, how='all').T
        index_value = results.index.get_level_values(1) if isinstance(results.index, pd.MultiIndex) else results.index
        results = results[index_value.notna()]
        results_dic[data] = results

    results_df = pd.concat(results_dic, axis=0, sort=False)
    results_df = results_df.T.drop_duplicates(keep="last").T
    return results_df

# If results are already generated, use the generated results:
if os.path.exists(metric_res_path / f"{metric}.tsv"):
    results = pd.read_table(metric_res_path / f"{metric}.tsv", index_col=0).T
else:
    results = get_method_metric(data_list=data_list,
                                methods=allMethods,
                                metric = metric)
    results = results.droplevel(-1)

#%% Heatmap
index_value = results.index.get_level_values(1) if isinstance(results.index, pd.MultiIndex) else results.index
results = results[index_value.notna()].astype(float)
# method must have over 5 non-na value
results = results.dropna(axis=1, thresh=4)
# results_NonNA = results.fillna(0.0)

col_anno = pd.DataFrame({"method": results.columns.str.split('_').str[0], 
                         # "mean": results.label
                         }, 
                         index=results.columns)    # col_dic = dict(zip(colored_col.unique(), wes_palette[:len(colored_col.unique())]))

row_anno = pd.read_table("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Cluster_Benchmark/Visualization/011_metric_results/samples_meta.tsv", 
                            index_col=0)

# Group the methods by SAMPLE (level0) or DOMAIN (level1), row_split changes correspondingly
col_anno["data_var"] = results.groupby(by=row_anno.dataset).mean().var().T
row_split = row_anno.dataset

box = anno_boxplot(results, cmap="cividis", height=30, legend=False)
bar = anno_barplot(results.mean(axis=1).apply(lambda x:round(x,2)), cmap="cividis", height=15, legend=False)
box.ylim = box.ylim[::-1]
row_ha = HeatmapAnnotation(
    mean=bar, #anno_boxplot(results, cmap="cividis", height=20, legend=False),
    tech=anno_simple(row_anno.technology, 
                           colors={"Visium": "#D8A499",  "STARmap+": "#5FBFF9", 
                                   "MERFISH": "#A4243B", "Xenium": "#E7BB41"}),
    data=anno_simple(row_anno.dataset, legend=False),
    dataset=anno_label(row_split, merge=True, rotation=30), 
    axis=0, 
    wgap = 0.2, 
    label_side="bottom", 
    # label_kws={'rotation':0,'horizontalalignment':'left','verticalalignment':'bottom'},
    orientation="right"
    )

"""
if isinstance(results.index, pd.MultiIndex):
    row_anno = pd.DataFrame({
        "data" : results.index.get_level_values(0),
        "samples" : results.index.get_level_values(1),
    }, index= results.index)

    # Group the methods by SAMPLE (level0) or DOMAIN (level1), row_split changes correspondingly
    col_boxplot = results.groupby(level=0).mean().T
    row_split = row_anno.data

    box = anno_boxplot(results, cmap="cividis", height=30, legend=False)
    box.ylim = box.ylim[::-1]
    row_ha = HeatmapAnnotation(sample=anno_label(row_split, merge=True, rotation=0), 
                               data=anno_simple(row_anno.data), 
                               # samples=anno_simple(row_anno.samples),
                               all_method_result=box, #anno_boxplot(results, cmap="cividis", height=20, legend=False),
                               axis=0, orientation='left')
"""

plt.figure(figsize=(13, 11))
col_ha = HeatmapAnnotation(method=anno_label(col_anno.method, merge=True,rotation=30),
                            methods=anno_simple(col_anno.method, legend=False),
                            data_var=anno_barplot(df=col_anno[["data_var"]], 
                                                       cmap = "plasma", height=15, legend=False),
                            axis=1,hgap=0.3, label_side="right")

cm = ClusterMapPlotter(data=results, 
top_annotation=col_ha, 
right_annotation=row_ha,
col_split = col_anno.method, 
col_split_order=col_anno.method.unique().tolist(),
col_cluster = True, 
row_split=row_split,
row_split_order=row_split.unique().tolist(), 
row_cluster=True, 
col_split_gap = 1, 
row_split_gap = 1,
show_rownames=False, 
legend_side = "right",
label = metric,
# xticklabels = True,
cmap = "viridis",
plot = True)

plt.suptitle(f"{metric}", y = 0.98)
plt.savefig(save_dir / f"{metric}_Hmp.png", bbox_inches="tight")

#%% Boxplots
def plot_method_boxplot(result, methods, metric, save_dir = None, nrow = None, ncol = None):
    import re
    import textwrap
    import math

    wes_palette = [
    # "#D8A499",  "#5FBFF9", "#A4243B", "#E7BB41",
    "#87B2C7", "#9B5094", "#6A994E", "#F1BB7B",
    "#FD6467", "#DF9D5A", "#F2D3BC", "#DAD6D6",
    "#2F343B", "#FAD77B", "#7D4E57", "#3E442B",
    ]

    newCol = []
    for s in result.columns:
        match = re.match(r'label(.*?) / config_(\d+)', s)
        if match:
            method = match.group(1)  # Extract the method name
            config = match.group(2)  # Extract the config number
            newCol.append(f'{method}_{config}')
        else:
            newCol.append(s)
    result.columns = newCol

    nrow = math.ceil(math.sqrt(len(methods))) if nrow is None else nrow
    ncol = math.ceil(len(methods) / nrow) if ncol is None else ncol
    fig, axs = plt.subplots(nrows=nrow, ncols=ncol, figsize=(7*ncol, 5*nrow))
    axf = axs.flatten()
    
    for i, method in enumerate(methods):
        data = result.filter(regex=f"^{method}*", axis=1)
        id_var = ["data", "sample"]
        data = pd.melt(data.reset_index(names=id_var), id_vars=id_var, var_name='method_instance', value_name=metric)
        data = data.dropna(axis=0)
        data["config"] = data['method_instance'].str.split('_').str[1]

        # plt.figure(figsize=(10, 4))  # Adjust figure size if needed
        sns.boxplot(data=data,  x = "data", y = metric, hue = "config", ax=axf[i],
                         dodge=True, gap=0.1, palette = wes_palette)
        sns.stripplot(data=data, x = "data", y = metric, ax=axf[i], hue = "config",
                    jitter=True, dodge=True, palette = wes_palette, legend=False, 
                    size=5, alpha = 0.8)
        handles, labels = axf[i].get_legend_handles_labels()
        axf[i].legend(handles[:], labels[:], title='configs', bbox_to_anchor=(1, 1.02), loc='upper left')
        axf[i].set_title(f"{method}: {metric}")

        labels = [textwrap.fill(label.get_text(), width=10) for label in axf[i].get_xticklabels()]
        axf[i].set_xticklabels(labels)
#        plt.gca().tick_params(axis='x', which='major', pad=5, length = 10)
        # plt.show()
        if save_dir is not None:
            axf[i].figure.savefig(save_dir / f"{method}_{metric}.png", bbox_inches="tight")
    
    fig.suptitle(f"{metric} for methods")
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig.show()

plot_method_boxplot(result=results, 
                    methods=["CellCharter", "STAGATE", "spaGCN", 
                             "SpaceFlow", "BANKSY", "SEDR", 
                             "seurat", "scanpy", "bass", 
                             "DRSC", "precast", "GraphST"], 
                    metric=metric,
                    save_dir=save_dir,
                    nrow=3)

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


# %%
