#%% Load packages
import os
from pathlib import Path
import sklearn.metrics
import pandas as pd
import re

import warnings
warnings.filterwarnings("ignore")

#%% Define methods and datasets

metric_list = [
    "adjusted_rand_score",
    "normalized_mutual_info_score",
    "fowlkes_mallows_score",
    "completeness_score",
    "homogeneity_score",
    "v_measure_score"
]

# metric = "v_measure_score"

data_list = [
    "libd_dlpfc",
    "visium_breast_cancer_SEDR",
    "visium_chicken_heart",
    "STARmap_plus",
    "abc_atlas_wmb_thalamus",
    "xenium-mouse-brain-SergioSalas",
    # "stereoseq_developing_Drosophila_embryos_larvae"
]

domain_path = Path("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Cluster_Benchmark/Visualization/010_results") 
result_path = Path("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Cluster_Benchmark/Visualization/011_metric_results") 
result_path.mkdir(parents=True, exist_ok=True)
dataset_path = Path("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Cluster_Benchmark/data")

#%% Define function to get metric results
def get_metric_results(metric, data_list, domain_path, result_path, dataset_path, metric_kw = {}, save=False):
    metric_res_dic = {}
    for data in data_list:
        data_path = domain_path / data

        # Get all observation.tsv from the results path
        if not any(file.endswith('.tsv') for file in os.listdir(data_path)):
            raise Exception("Observation tsv not saved yet! Run 010_Save_Results.py to generate the results first!")
        samples = [m for m in os.listdir(dataset_path / data) if os.path.isdir(os.path.join(dataset_path, data, m)) and not m.startswith(".")]
        obss = [pd.read_table(data_path / f"obs_{sample}.tsv", index_col=0) for sample in samples]

        # Get metric results from each sample and columns
        for i, obs in enumerate(obss):
            # Leave only columns indicating clusteirng results
            domain = obs.filter(regex="^label*", axis=1)
            newCol = []

            metric_res_dic[samples[i]] = pd.Series([getattr(sklearn.metrics, metric)(
                    domain["label"].astype("category").cat.codes, 
                    domain[col].astype("category").cat.codes, 
                    **metric_kw
                    ) 
                    for col in domain.columns[domain.columns!="label"]], 
                    index=domain.columns[domain.columns!="label"])
    metric_df = pd.concat(metric_res_dic, axis=1, sort=False)

    # Rename columns by deleting the "/config_" in the name

    for s in metric_df.index:
        if s.split("label")[1] != "":
            newCol.append(s.split("label")[1])
        else:
            newCol.append(s)
    metric_df.index = newCol
    
    if save:
        metric_df.to_csv(result_path / f"{metric}.tsv", sep="\t", index_label="")

    return metric_df


#%% Generating all the metric results
for metric in metric_list:
    met_df_list = get_metric_results(metric=metric, 
                                     data_list=data_list, 
                                     dataset_path=dataset_path,
                                     result_path=result_path,
                                     domain_path=domain_path,
                                     save=True)

# %% Save samples dataframe
meta_list = []
for data in data_list:
    import json

    data_path = dataset_path / data
    samples = [m for m in os.listdir(dataset_path / data) if os.path.isdir(os.path.join(dataset_path, data, m)) and not m.startswith(".")]
    with open(data_path / "experiment.json", "r") as file:
        info = json.load(file)

    sample_df = pd.read_table(data_path / "samples.tsv")[["n_clusters"]]
    sample_df.index = samples
    sample_df["dataset"] = data
    sample_df["technology"] = info["technology"]
    sample_df["n_genes"] = [len(pd.read_table(data_path / sample / "features.tsv")) for sample in samples]
    sample_df["n_cells"] = [len(pd.read_table(data_path / sample / "observations.tsv")) for sample in samples]
    meta_list.append(sample_df)

meta_df = pd.concat(meta_list, axis=0, sort=False)
meta_df.to_csv(result_path / "samples_meta.tsv", sep="\t", index_label="")
# %%
