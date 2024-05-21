#%% Load packages
import os
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
from PyComplexHeatmap import *

import sklearn.metrics
import pandas as pd
import itertools
import re

import warnings
warnings.filterwarnings("ignore")

#%% Define methods and datasets

metric = "adjusted_rand_score"

# data = "libd_dlpfc"
# data = "visium_breast_cancer_SEDR"
# data = "visium_chicken_heart"
data = "STARmap_plus"
# data = "abc_atlas_wmb_thalamus"
# data = "xenium-mouse-brain-SergioSalas"
# data = "stereoseq_developing_Drosophila_embryos_larvae"

domain_path = Path("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Cluster_Benchmark/Visualization/010_results") / data
result_path = Path("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Cluster_Benchmark/Visualization/012_Domain_Cross") / data

dataset_path = Path("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Cluster_Benchmark/data")
data_path = dataset_path / data
samples = [m for m in os.listdir(data_path) if os.path.isdir(os.path.join(data_path, m)) and not m.startswith(".")]

#%% Get adata object
if not any(file.endswith('.tsv') for file in os.listdir(domain_path)):
      raise Exception("Observation tsv not saved yet! Run 010_Save_Results.py to generate the results first!")
obss = [pd.read_table(domain_path / f"obs_{sample}.tsv", index_col=0) for sample in samples]

res_cross_list = []
for i, obs in enumerate(obss):
    domain = obs.filter(regex="^label*", axis=1)
    newCol = []
    for s in domain.columns:
        if s.split("label")[1] != "":
            newCol.append(s.split("label")[1])
        else:
            newCol.append(s)
    domain.columns = newCol
    res_comb = itertools.combinations(domain.columns, 2)
    res_cross = pd.DataFrame(columns=domain.columns, index=domain.columns, data=0)

    for col1, col2 in res_comb:
        ari_score = sklearn.metrics.adjusted_rand_score(domain[col1].astype("category").cat.codes, 
                                        domain[col2].astype("category").cat.codes)
        res_cross.loc[col1, col2] = res_cross.loc[col2, col1] = ari_score

    save_path = result_path / samples[i]
    save_path.mkdir(exist_ok=True, parents=True)
    res_cross.to_csv(save_path / f"{metric}.tsv", sep="\t", index_label="")
    res_cross_list.append(res_cross)


#%% plot and save spatial data
wes_palette = [
    "#D8A499",  "#5FBFF9", "#A4243B", "#E7BB41",
    "#2F343B", "#FAD77B", "#7D4E57", "#3E442B",
    "#87B2C7", "#9B5094", "#6A994E", "#F1BB7B",
    "#FD6467", "#DF9D5A", "#F2D3BC", "#DAD6D6",
]
for i, result in enumerate(res_cross_list):

    result.fillna(0.0, inplace=True)
    index_value = result.index.get_level_values(1) if isinstance(result.index, pd.MultiIndex) else result.index
    result = result[index_value.notna()].astype(float)

    method_anno = pd.DataFrame({"method": result.columns.str.split('_').str[0], 
                            "label": result.label}, 
                            index=result.index)    # col_dic = dict(zip(colored_col.unique(), wes_palette[:len(colored_col.unique())]))

    plt.figure(figsize=(9, 8))
    col_ha = HeatmapAnnotation(label=anno_label(method_anno.method, merge=True,rotation=30),
                                method=anno_simple(method_anno.method, legend=False),
                                #GT_Metric=anno_barplot(method_anno.label, cmap = "plasma", height=30),
                                axis=1,)
    cm = ClusterMapPlotter(data=result, 
    top_annotation=col_ha,
    left_annotation=HeatmapAnnotation(label=anno_label(method_anno.method, merge=True), 
                                      method=anno_simple(method_anno.method, legend=False), axis=0), 
    col_split = method_anno.method, 
    col_cluster = False, 
    row_split=method_anno.method, 
    row_cluster=False, 
    col_split_gap = 1, 
    show_rownames=False, 
    show_colnames=False,
    label = metric,
    cmap = "PuBu")

    plt.suptitle(f"{metric}: {data}_{samples[i]}", y = 0.98, fontsize = 16)
    save_path = result_path / samples[i]
    save_path.mkdir(exist_ok=True, parents=True)
    plt.savefig(save_path / f"{metric}_cross_Hmp.png", bbox_inches="tight")
# %% Plot a specific one 
result = pd.read_table("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Cluster_Benchmark/Visualization/012_Domain_Cross/STARmap_plus/well07/adjusted_rand_score.tsv",
                       index_col=0)

wes_palette = [
    "#D8A499",  "#5FBFF9", "#A4243B", "#E7BB41",
    "#2F343B", "#FAD77B", "#7D4E57", "#3E442B",
    "#87B2C7", "#9B5094", "#6A994E", "#F1BB7B",
    "#FD6467", "#DF9D5A", "#F2D3BC", "#DAD6D6",
]

result.fillna(0.0, inplace=True)
index_value = result.index.get_level_values(1) if isinstance(result.index, pd.MultiIndex) else result.index
result = result[index_value.notna()].astype(float)

method_anno = pd.DataFrame({"method": result.columns.str.split('_').str[0], 
                        "label": result.label}, 
                        index=result.index)    # col_dic = dict(zip(colored_col.unique(), wes_palette[:len(colored_col.unique())]))

plt.figure(figsize=(9, 8))
col_ha = HeatmapAnnotation(label=anno_label(method_anno.method, merge=True,rotation=30),
                            method=anno_simple(method_anno.method, legend=False),
                            #GT_Metric=anno_barplot(method_anno.label, cmap = "plasma", height=30),
                            axis=1,)
cm = ClusterMapPlotter(data=result, 
top_annotation=col_ha,
left_annotation=HeatmapAnnotation(label=anno_label(method_anno.method, merge=True), 
                                    method=anno_simple(method_anno.method, legend=False), axis=0), 
col_split = method_anno.method, 
col_cluster = False, 
row_split=method_anno.method, 
row_cluster=False, 
col_split_gap = 1, 
show_rownames=False, 
# row_names_side="left",
label = 'ARI',
xticklabels = True,
cmap = "PuBu")

# %%
