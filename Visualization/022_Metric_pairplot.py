#%% Define packages
import os
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from functools import reduce


#%% Define get ARI function
metrics = [
    "adjusted_rand_score",
    "NMI",
    "fowlkes_mallows_score",
    "homogeneity_score",
    "completeness_score",
    "v_measure_score"
]

metric_names = ["ARI", "NMI", "FMI", "HOM", "COM", "VMS"]
metric_res_path = Path("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Cluster_Benchmark/Visualization/011_metric_results")

res_list = []

# Collect all metric values
for i, metric in enumerate(metrics):
    # If results are already generated, use the generated results:
    if os.path.exists(metric_res_path / f"{metric}.tsv"):
        results = pd.read_table(metric_res_path / f"{metric}.tsv", index_col=0).T
    else:
        raise ValueError("No metric results found! Please generate the metric results first")

    results = results.reset_index(names="sample")
    results = pd.melt(results, id_vars = "sample", value_name=metric_names[i], var_name="methods")
    res_list.append(results)

all_metrics = reduce(lambda left, right: pd.merge(left, right, on=['sample', 'methods'], how="outer"), res_list)

sample_meta = pd.read_table("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Cluster_Benchmark/Visualization/011_metric_results/samples_meta.tsv", 
                            index_col=0)["dataset"]
all_metrics["data"] = sample_meta[all_metrics["sample"]].to_list()
all_metrics["methods"] = all_metrics['methods'].str.split('_').str[0]
all_metrics = all_metrics.dropna(subset=metric_names, how="all")

#%% Pairplot
pg1 = sns.pairplot(all_metrics, hue="data", corner=True)
pg2 = sns.pairplot(all_metrics, hue="methods", corner=False)

pg1.savefig("data_pair.png", transparent=True)
pg2.savefig("method_pair.png", transparent=True)

# %%
