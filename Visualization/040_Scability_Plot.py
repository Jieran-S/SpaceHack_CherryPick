#%% Load packages
import os
import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
from PyComplexHeatmap import *

import warnings
warnings.filterwarnings("ignore")

#%% define function to get results
def get_data_usability(methods, data_path, 
                       samples, 
                       object: str | None = "memory"):
    import os
    import pandas as pd

    match object:
        case "memory":
            object = "max_uss"
        case "runtime":
            object = "s"

    def extract_meta(file_path, object):
        meta_res = pd.read_table(file_path)
        return meta_res.loc[0, object]

    results_dict = {}
    for sample in samples:
        sample_result = {}
        avail_methods = set(methods).intersection(os.listdir(os.path.join(data_path, sample)))
        if len(avail_methods) == 0:
            continue
        else:
            for method in avail_methods:
                method_dir = os.path.join(data_path, sample, method)
                if "benchmark_method.txt" in os.listdir(method_dir):
                    sample_result[method] = extract_meta(os.path.join(method_dir, "benchmark_method.txt"), 
                                                 object)
                else:
                    configs = [x for x in os.listdir(method_dir) if x.startswith("config")]
                    for config in configs:
                        if "benchmark_method.txt" in os.listdir(os.path.join(method_dir, config)):
                            config_name = config.split("config_")[1]
                            if "leiden" in config:
                                mth_ins = f"{method}-leiden_{config_name}"
                            elif "louvain" in config:
                                mth_ins = f"{method}-louvain_{config_name}"
                            else:
                                mth_ins = f"{method}_{config_name}"
                            sample_result[mth_ins] = extract_meta(os.path.join(method_dir, config, "benchmark_method.txt"), 
                                                          object)
        results_dict[sample] = pd.DataFrame.from_dict(sample_result, orient="index")
    meta_df = pd.concat(results_dict, axis=1)
    meta_df = meta_df.sort_index(axis=0)
    return meta_df

def get_method_usability(dataset_path, data_list, methods, object=None, prior_dir=None):
    from pathlib import Path
    import pandas as pd
    import os

    if prior_dir is None:
        prior_dir = Path(f"/work/PRTNR/CHUV/DIR/rgottar1/spatial/Cluster_Benchmark/Visualization/040_Scability")

    results_dic = {}
    for data in data_list:
#        if os.path.isfile(prior_dir / f"{metric}_{data}.tsv"):
#            results = pd.read_table(prior_dir / f"{metric}_{data}.tsv", index_col=0)
#            results_list.append(results)
#        else: 
        data_path = dataset_path / data
        samples = [m for m in os.listdir(data_path) if os.path.isdir(os.path.join(data_path, m)) and not m.startswith(".")]
        results = get_data_usability(methods= methods,
                        data_path=data_path,
                        samples=samples,
                        object=object)

        # if results.shape[1] > 1:
        results = results.dropna(axis=1, how='all').T
        index_value = results.index.get_level_values(1) if isinstance(results.index, pd.MultiIndex) else results.index
        results = results[index_value.notna()]
        results_dic[data] = results

    results_df = pd.concat(results_dic, axis=0, sort=False)
    results_df = results_df.T.drop_duplicates(keep="last").T
    results_df = results_df.droplevel(-1)
    return results_df

#%% Define methods and datasets
wes_palette_final = [
  "#710421",  # dark Red
  "#Ff5100",  # International Orange
  "#836247",  # coffee light
  "#F9C846",  # Bright Mustard Yellow
  "#34568B",  # Dark Blue
  "#175620",  # Parsley green
  "#8e67d8",  # Medium Purple
  "#1d8af9",  # Dodger Blue
  "#DC143C",  # Crimson
  "#005A51",  # Deep Teal
  "#381756",  # grape Purple
  "#9c9353",  # Corn dark yellow
  "#99D6EA",  # Brighter Sky Blue
  "#293037",  # Space gray
  "#69e256",  # Pastel Green
  "#C44E55",  # Bright Deep Red
]

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

dataset_path = Path("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Cluster_Benchmark/data/")
save_dir = Path(f"/work/PRTNR/CHUV/DIR/rgottar1/spatial/Cluster_Benchmark/Visualization/040_Scability")
data_list = [
    "libd_dlpfc",
    "visium_breast_cancer_SEDR",
    "visium_chicken_heart",
    "STARmap_plus",
    "abc_atlas_wmb_thalamus",
    "xenium-mouse-brain-SergioSalas",
    # "stereoseq_developing_Drosophila_embryos_larvae"
]
tech_dic = {
    "libd_dlpfc":"Visium",
    "visium_breast_cancer_SEDR":"Visium",
    "visium_chicken_heart":"Visium",
    "STARmap_plus":"Starmap+",
    "abc_atlas_wmb_thalamus":"MERFISH",
    "xenium-mouse-brain-SergioSalas":"Xenium",
    "stereoseq_developing_Drosophila_embryos_larvae":"Stereoseq"
}

#%% Get adata object
if os.path.exists("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Cluster_Benchmark/Visualization/040_Usability/scalability.tsv"):
    data = pd.read_table("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Cluster_Benchmark/Visualization/040_Usability/scalability.tsv",
                            index_col=0)
else:
    memory_df = get_method_usability(dataset_path=dataset_path,
                                data_list=data_list,
                                methods=allMethods,
                                object="memory")

    runtime_df = get_method_usability(dataset_path=dataset_path,
                                data_list=data_list,
                                methods=allMethods,
                                object="runtime")

#%% save results
memory_df.to_csv(save_dir / f"memory.tsv", sep="\t", index_label=False)
runtime_df.to_csv(save_dir / f"runtime.tsv", sep="\t", index_label=False)

#%% plot usability visualization: Overview heatmap
results = runtime_df
# object = "memory"
object = "runtime"

index_value = results.index.get_level_values(1) if isinstance(results.index, pd.MultiIndex) else results.index
results = results[index_value.notna()].astype(float)
# results_NonNA = results.fillna(0.0)

col_anno = pd.DataFrame({"method": results.columns.str.split('/').str[0], 
                         # "mean": results.label
                         }, 
                         index=results.columns)    # col_dic = dict(zip(colored_col.unique(), wes_palette[:len(colored_col.unique())]))

row_anno = pd.DataFrame({
    "data" : results.index.get_level_values(0),
    "samples" : results.index.get_level_values(1),
}, index= results.index)

# Group the methods by SAMPLE (level0) or DOMAIN (level1), row_split changes correspondingly
col_boxplot = results.groupby(level=0).mean().T
row_split = row_anno.data

bar = anno_barplot(results.mean(axis=1), cmap="cividis", height=30, legend=False)

row_ha = HeatmapAnnotation(sample=anno_label(row_split, merge=True, rotation=60), 
                            data=anno_simple(row_anno.data, legend=False), 
                            # samples=anno_simple(row_anno.samples),
                            # average=bar,
                            axis=0, orientation='left')

plt.figure(figsize=(13, 11))
col_ha = HeatmapAnnotation(method=anno_label(col_anno.method, merge=True,rotation=30),
                            methods=anno_simple(col_anno.method),
                            # Domain_metric=anno_boxplot(col_boxplot, cmap = "plasma", height=30, legend=False),
                            mean=anno_barplot(results.mean(axis=0), cmap="plasma", legend=False),
                            axis=1,)
cm = ClusterMapPlotter(data=results, 
top_annotation=col_ha, 
left_annotation=row_ha,
col_split = col_anno.method, 
col_cluster = False, 
row_split=row_split, 
row_cluster=False, 
col_split_gap = 1, 
show_rownames=False, 
# row_names_side="left",
label = object,
# xticklabels = True,
# cmap = "Wistia"
)

plt.suptitle(f"{object}", y = 0.98)
plt.savefig(save_dir / f"{object}_Hmp.png", bbox_inches="tight")

#%% data formatting
load_exist = True
if "scalability.tsv" in os.listdir(save_dir) and load_exist:
    data = pd.read_table(save_dir / "scalability.tsv", index_col=0)
else:
    id_var = ["data", "sample"]
    data = pd.melt(memory_df.reset_index(names=id_var), id_vars=id_var, 
                var_name='method_instance', value_name="memory_mb")
    runtime_data = pd.melt(runtime_df.reset_index(names=id_var), id_vars=id_var, 
                var_name='method_instance', value_name="runtime_s")
    data = data.join(runtime_data["runtime_s"])
    data["method"] = data['method_instance'].str.split('_').str[0]

    n_cell = []
    n_gene = []

    for i, row in data.iterrows():
        sample_dir = dataset_path / row["data"] / row["sample"]
        n_cell.append(len(pd.read_table(sample_dir / "observations.tsv")))
        n_gene.append(len(pd.read_table(sample_dir / "features.tsv")))

    data["n_cell"] = n_cell
    data["n_gene"] = n_gene
    data["technology"] = data["data"].map(tech_dic)

    metric = "ARI"
    metric_res = pd.read_table(f"/work/PRTNR/CHUV/DIR/rgottar1/spatial/Cluster_Benchmark/Visualization/021_CrossSample/{metric}.tsv",
                            index_col=[0,1])
    metric_res = pd.melt(metric_res.reset_index(names=id_var), id_vars=id_var, 
                    var_name='method_instance', value_name=metric)

    data["runtime(min)"] = data["runtime_s"]/60
    data["memory(GB)"] = data["memory_mb"]/1000

    data = data.join(metric_res[metric])
    data = data.dropna(axis=0)

    data.to_csv(save_dir / "scalability.tsv", sep= "\t", index_label=False)

#%% Plotting scatterplot with density
data["runtime(min)"] = data["runtime_s"]/60
data["memory(GB)"] = data["memory_mb"]/1000

color_dict = {}
for i, method in enumerate(data.method.unique()):
    color_dict[method] = wes_palette_final[i]

plt.figure(figsize=(7, 7))  # Adjust figure size if needed
ax = sns.scatterplot(data=data, x="runtime(min)", y = "memory(GB)", color="k", 
                     alpha = 1, linewidth=0, style="data", 
                     hue = "method", palette=color_dict, s=60)
sns.kdeplot(data=data, x="runtime(min)", y = "memory(GB)", ax= ax, alpha = 0.3)
ax.axvline(x=1440, color='b', linestyle='--')
ax.axhline(y=500, color='b', linestyle='--')
ax.set_xscale("log", base=2)
ax.set_yscale("log", base=2)
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[:], labels[:], title='method', bbox_to_anchor=(1, 1.02), loc='upper left')

# ax.legend(handles[:], labels[:], title='method', bbox_to_anchor=(1, 1.02), loc='upper left')
plt.title(f"Usability")
# plt.xticks(rotation=90)  # Rotate x-axis labels for better readability
# plt.show()
plt.savefig(save_dir / f"Runtime_VS_Memory.png", bbox_inches="tight")

#%% Sample-size based runtime/memory usage
data["runtime(min)"] = data["runtime_s"]/60
data["memory(GB)"] = data["memory_mb"]/1000

fig, axs = plt.subplots(nrows=4, ncols=1, figsize=(14*1, 5*4))
axf = axs.flatten()

data_list = [
    data.loc[data["technology"] == "Visium", :],
    data.loc[data["technology"] == "Visium", :],
    data.loc[data["technology"] != "Visium", :],
    data.loc[data["technology"] != "Visium", :],
]

for i, data_sub in enumerate(data_list):
    y_lab = ["memory(GB)", "runtime(min)", "memory(GB)", "runtime(min)"]
    sns.lineplot(data=data_sub, x="n_cell", y = y_lab[i], color="k", 
                    alpha = 1, 
                    ax=axf[i], marker="o",
                    hue = "method", palette=wes_palette_final)

    for line, name in zip(axf[i].lines, data_sub["method"].unique()):
        y = line.get_ydata()[-1]
        x = line.get_xdata()[-1]
        if not np.isfinite(y):
            y=next(reversed(line.get_ydata()[~line.get_ydata().mask]),float("nan"))
        if not np.isfinite(y) or not np.isfinite(x):
            continue     
        text = axf[i].annotate(name,
                xy=(x*1.01, y),
                xytext=(0, 0),
                color=line.get_color(),
                xycoords=(axf[i].get_xaxis_transform(),
                    axf[i].get_yaxis_transform()),
                textcoords="offset points")

    handles, labels = axf[i].get_legend_handles_labels()
    axf[i].legend(handles[:], labels[:], title='method', bbox_to_anchor=(1, 1.02), loc='upper left')
    axf[i].set_title(["Visium: Memory", "Visium: Runtime", 
                      "Non-Visium: Memory", "Non-Visium: Runtime"][i])

fig.suptitle(f"Runtime and Memory in Visium and NonVisium")
fig.tight_layout(rect=[0, 0.03, 1, 0.95])
fig.show()

# %% Runtime and accuracy -> Per-method
fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(10*2, 5*2))
axf = axs.flatten()

data_list = [
    data.loc[data["technology"] == "Visium", :],
    data.loc[data["technology"] == "Visium", :],
    data.loc[data["technology"] != "Visium", :],
    data.loc[data["technology"] != "Visium", :],
]

for i, data_sub in enumerate(data_list):
    x_lab = ["memory(GB)", "runtime(min)", "memory(GB)", "runtime(min)"]
    for j, met in enumerate(data_sub.method.unique()):
        sns.regplot(data=data_sub[data_sub.method==met], 
                    x=x_lab[i], 
                    y = metric, 
                    label = met, 
                    ax = axf[i],
                    ci=None,
                    color=wes_palette_final[j])

    handles, labels = axf[i].get_legend_handles_labels()
    axf[i].legend(handles[:], labels[:], title='method', bbox_to_anchor=(1, 1.02), loc='upper left')
    axf[i].set_title(["Visium: Memory", "Visium: Runtime", 
                      "Non-Visium: Memory", "Non-Visium: Runtime"][i])

fig.suptitle(f"ARI vs Memory/Runtime")
fig.tight_layout(rect=[0, 0.03, 1, 0.98])
fig.show()

#%% Product of runtime and memory
data["runtime_x_memory"] = data["runtime(min)"] * data["memory(GB)"]
data["n_gene_x_n_cell"] = data["n_gene"] * data["n_cell"]


data_list = [
    data.loc[data["technology"] != "Visium", :],
    data.loc[data["technology"] == "Visium", :],
    # data.loc[data["technology"] == "Visium", :],
    # data.loc[data["technology"] != "Visium", :],
]

color_dict = {}
for i, method in enumerate(data.method.unique()):
    color_dict[method] = wes_palette_final[i]

#%% accuracy relationship 
nrow = 2
ncol = 1
fig, axs = plt.subplots(nrows=nrow, ncols=ncol, figsize=(10*ncol, 5*nrow))
axf = axs.flatten()

for i, data_sub in enumerate(data_list):
    #y_lab = ["memory(GB)", "runtime(min)", "memory(GB)", "runtime(min)"]
    for j, met in enumerate(data_sub.method.unique()):
        sns.regplot(data=data_sub[data_sub.method==met], 
                    x="runtime_x_memory", 
                    y = "ARI", 
                    label = met, 
                    ax = axf[i],
                    ci=None,
                    robust=True,
                    scatter_kws={"alpha":0.3, "s":5}, 
                    color=color_dict[met])
    axf[i].set_xscale("log")
    axf[i].set_xlabel("Runtime x Memory")

    for line, name in zip(axf[i].lines, data_sub["method"].unique()):
        y = line.get_ydata()[-1]
        x = line.get_xdata()[-1]
        if not np.isfinite(y):
            y=next(reversed(line.get_ydata()[~line.get_ydata().mask]),float("nan"))
        if not np.isfinite(y) or not np.isfinite(x):
            continue     
        text = axf[i].annotate(name,
                xy=(x*1.01, y),
                xytext=(0, 0),
                color=line.get_color(),
                xycoords=(axf[i].get_xaxis_transform(),
                    axf[i].get_yaxis_transform()),
                textcoords="offset points")

    handles, labels = axf[i].get_legend_handles_labels()
    if i==4:
        axf[i].legend(handles[:], labels[:], title='method', bbox_to_anchor=(1, 1.02), loc='upper left')
    axf[i].set_title(["Visium", "Non-Visium"][i])

#fig.suptitle(f"Runtime and Memory in Visium and NonVisium")
#fig.tight_layout(rect=[0, 0.03, 1, 0.95])
fig.show()

# %% Cell and Gene relationship 
fig, axs = plt.subplots(nrows=4, ncols=1, figsize=(7*1, 5*4))
axf = axs.flatten()

y_lab = ["memory(GB)", "runtime(min)", "memory(GB)", "runtime(min)"]

for i, data_sub in enumerate(data_list):
    sns.lineplot(data=data_sub, x="n_gene_x_n_cell", 
                 y = y_lab[i], 
                 color="k", 
                    alpha = 1, 
                    ax=axf[i], marker="o",
                    hue = "method",
                    palette=color_dict)

    for line, name in zip(axf[i].lines, data_sub["method"].unique()):
        y = line.get_ydata()[-1]
        x = line.get_xdata()[-1]
        if not np.isfinite(y):
            y=next(reversed(line.get_ydata()[~line.get_ydata().mask]),float("nan"))
        if not np.isfinite(y) or not np.isfinite(x):
            continue     
        text = axf[i].annotate(name,
                xy=(x*1.01, y),
                xytext=(0, 0),
                color=line.get_color(),
                xycoords=(axf[i].get_xaxis_transform(),
                    axf[i].get_yaxis_transform()),
                textcoords="offset points")

    axf[i].set_xscale("log")
    handles, labels = axf[i].get_legend_handles_labels()
    axf[i].legend(handles[:], labels[:], title='method', bbox_to_anchor=(1, 1.02), loc='upper left')
    axf[i].set_title(["Visium: Memory", "Visium: Runtime", 
                      "Non-Visium: Memory", "Non-Visium: Runtime"][i])

fig.suptitle(f"Runtime and Memory in Visium and NonVisium")
fig.tight_layout(rect=[0, 0.03, 1, 0.95])
fig.show()
# %%
