#%% Load packages
import os
from pathlib import Path
import matplotlib.pyplot as plt
import scanpy as sc

import warnings
warnings.filterwarnings("ignore")

#%% define function to get results
def get_anndata(folder, qc=False):
    # Untested template
    import anndata as ad
    import pandas as pd
    import scipy as sp
    from pathlib import Path

    folder = Path(folder) / "qc" if qc else Path(folder)

    X = sp.io.mmread(folder / "counts.mtx")
    if sp.sparse.issparse(X):
        X = X.tocsr()
    observations = pd.read_table(folder / "observations.tsv", index_col=0)
    features = pd.read_table(folder / "features.tsv", index_col=0)
    coordinates = (
        pd.read_table(folder / "coordinates.tsv", index_col=0)
        .loc[observations.index, :]
        .to_numpy()
    )

    adata = ad.AnnData(
        X=X, obs=observations, var=features, obsm={"spatial": coordinates}
    )

    return adata

def get_domains(data_path, methods):
    import pandas as pd
    import os 
    from pathlib import Path

    res_list = []
    # print(samples)
    for sample in samples:
        adata = get_anndata(Path(data_path) / sample, qc=True)
        adata.uns["sample"] = sample
        res_df = adata.obs
        if os.path.exists(os.path.join(data_path, sample, "labels.tsv")):
            label_df = pd.read_table(os.path.join(data_path, sample, "labels.tsv"), index_col=0)["label"].astype("category")
            res_df = res_df.join(label_df)

        avail_methods = set(methods).intersection(os.listdir(os.path.join(data_path, sample)))
        if len(avail_methods) != 0:
            for method in avail_methods:
                method_dir = os.path.join(data_path, sample, method)
                if "domains.tsv" in os.listdir(method_dir):
                    res_df_mth = pd.read_table(os.path.join(method_dir, "domains.tsv"), index_col=0)["label"].astype("category")
                    res_df = res_df.join(res_df_mth,  rsuffix=method)
                else:
                    configs = [x for x in os.listdir(method_dir) if x.startswith("config")]
                    for config in configs:
                        if "domains.tsv" in os.listdir(os.path.join(method_dir, config)):
                            ari_df_mth = pd.read_table(os.path.join(method_dir, config, "domains.tsv"), index_col=0)["label"].astype("category")
                            res_df = res_df.join(ari_df_mth,  rsuffix=f"{method} / {config}")
        res_df = res_df.T.drop_duplicates(keep="last").T
        cols = [lab for lab in res_df.columns if lab.startswith("label") and lab != "label"]
        res_df[cols] = res_df[cols].stack().astype("category").unstack()
        adata.obs = res_df
        res_list.append(adata)
    return res_list

#%% Define methods and datasets
result_path = Path("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Cluster_Benchmark/Visualization/010_results")

# data_name = "libd_dlpfc"
# data_name = "visium_breast_cancer_SEDR"
# data_name = "visium_chicken_heart"
# data_name = "STARmap_plus"
# data_name = "abc_atlas_wmb_thalamus"
data_name = "xenium-mouse-brain-SergioSalas"
# data_name = "stereoseq_developing_Drosophila_embryos_larvae"
data_path = result_path / data_name

#%% Get adata object

ad_list = [sc.read_h5ad(data_path / x) for x in os.listdir(data_path) if x.endswith("h5ad")]

#%% plot and save spatial data
for adata in ad_list:
        # adata.obs = adata.obs.T.drop_duplicates().T
        save_dir = Path(f"/work/PRTNR/CHUV/DIR/rgottar1/spatial/Cluster_Benchmark/Visualization/030_Spatialplots/{data_name}")
        save_dir.mkdir(exist_ok=True, parents=True)
        with plt.rc_context():  # Use this to set figure params like size and dpi
                sc.pl.spatial(adata, 
                        color=[lab for lab in adata.obs.columns if lab.startswith("label")], 
                        spot_size = 35, 
                        ncols = 7,
                        # save = str(save_dir / f"{adata.uns['sample']}_{data_name}.png"),
                        show = False)
                plt.savefig(save_dir / f"{adata.uns['sample']}_{data_name}.png", bbox_inches="tight")
                plt.show()

# %% all ground truth from all samples
data_list = [
     "libd_dlpfc",
     "visium_breast_cancer_SEDR",
     "visium_chicken_heart",
     "STARmap_plus",
     "abc_atlas_wmb_thalamus", 
     "xenium-mouse-brain-SergioSalas"
]

#%%
ad_old = sc.read_h5ad("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Cluster_Benchmark/Visualization/010_results_1204/visium_breast_cancer_SEDR/adata_visium_breast_cancer_SEDR.h5ad")

import sklearn.metrics
import pandas as pd
import re

obs = ad_old.obs
# metric = "adjusted_rand_score"
domain = obs.filter(regex="^label*", axis=1)
newCol = []

metric_res_dic = pd.Series([sklearn.metrics.v_measure_score(
        domain["label"].astype("category").cat.codes, 
        domain[col].astype("category").cat.codes) 
        for col in domain.columns[domain.columns!="label"]], 
        index=domain.columns[domain.columns!="label"])

old_ari = pd.DataFrame(metric_res_dic)
old_ari.columns = ["ARI_4"]
method = [s.split(" / ")[0] for s in old_ari.index]
method = [s.split("label")[1] for s in method]
old_ari["method"] = method

old_ari_mean = old_ari.groupby("method").mean()

new_ari = pd.read_table("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Cluster_Benchmark/Visualization/011_metric_results/v_measure_score.tsv", 
                        index_col=0)[["visium_breast_cancer_SEDR"]]
new_method = [s.split("_")[0] for s in new_ari.index]
new_ari["method"] = new_method
new_ari_mean = new_ari.groupby("method").mean()
new_ari_mean.columns = ["ARI_20"]

ari_SEDR = old_ari_mean.join(new_ari_mean, how="inner")
ari_SEDR["method"] = ari_SEDR.index 
ari_SEDR = ari_SEDR.melt(id_vars="method", value_name="V Measure")
ari_SEDR["Cluster Number"] = [s.split("ARI_")[1] for s in ari_SEDR.variable]

import seaborn as sns

fig, axf = plt.subplots(nrows=1, ncols=1, figsize=(10, 5))

sns.lineplot(data=ari_SEDR, x = "method", y = "V Measure", 
             hue = "Cluster Number", marker="o", sort= False,
             ax=axf)
axf.set_xticklabels(axf.get_xticklabels(), rotation=45)
handles, labels = axf.get_legend_handles_labels()
axf.legend(handles[:], labels[:], title='method', bbox_to_anchor=(0.5, 1.02), loc='upper center')

plt.show()

# %% Plot spatial plot with CHAOS score and ARI
import pandas as pd
adata = sc.read_h5ad("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Cluster_Benchmark/Visualization/010_results/xenium-mouse-brain-SergioSalas/adata_region_main.h5ad")
CHAOS_score = pd.read_table("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Cluster_Benchmark/Visualization/011_metric_results/CHAOS.tsv", index_col=0)["region_main"]
ARI_score = pd.read_table("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Cluster_Benchmark/Visualization/011_metric_results/adjusted_rand_score.tsv", index_col=0)["region_main"]

CHAOS_score["label"] = 0
ARI_score["label"] = 0

labels = [lab for lab in adata.obs.columns if lab.startswith("label")]
label_match =["label"] + [s.replace("label", "") for s in labels[1:]]

metric_df = pd.DataFrame({"CHAOS" : CHAOS_score,
                          "ARI": ARI_score,
                          "Columns": pd.Series(labels, index=label_match)})
metric_df = metric_df.sort_values(by = "ARI", ascending=False)

with plt.rc_context():  # Use this to set figure params like size and dpi
    sc.pl.spatial(adata, 
            color=metric_df["Columns"], 
            spot_size = 35, 
            ncols = 7,
            title = [f"{lab} CHAOS:{CHAOS_score[lab]} ARI:{ARI_score[lab]}" for lab in metric_df.index],
            # save = str(save_dir / f"{adata.uns['sample']}_{data_name}.png"),
            show = False)
    plt.show()

labs = ["label", "labelSEDR_2", "labelCellCharter_1"]
labs_mat = ["label", "SEDR_2", "CellCharter_1"]
sc.pl.spatial(adata, 
            color=labs, 
            spot_size = 35, 
            title = [f"{lab} CHAOS:{CHAOS_score[lab]} ARI:{ARI_score[lab]}" for lab in labs_mat],
            # save = str(save_dir / f"{adata.uns['sample']}_{data_name}.png"),
            show = False)
# %%
