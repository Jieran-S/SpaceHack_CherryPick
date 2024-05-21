#%% Load packages
import os
from pathlib import Path
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

def get_result_anndata(data_path, methods):
    import pandas as pd
    import os 
    from pathlib import Path

    res_list = []
    samples = [m for m in os.listdir(data_path) if os.path.isdir(os.path.join(data_path, m)) and not m.startswith(".")]
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
                            config_name = config.split("config_")[1]
                            if "leiden" in config:
                                suffix = f"{method}-leiden_{config_name}"
                            elif "louvain" in config:
                                suffix = f"{method}-louvain_{config_name}"
                            else:
                                suffix = f"{method}_{config_name}"
                            res_df = res_df.join(ari_df_mth,  rsuffix=suffix)
        # res_df = res_df.T.drop_duplicates(keep="last").T
        cols = [lab for lab in res_df.columns if lab.startswith("label") and lab != "label"]
        res_df[cols] = res_df[cols].stack().astype("category").unstack()
        adata.obs = res_df
        res_list.append(adata)
    return res_list

#%% Define methods and datasets
all_methods =["BayesSpace", 
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
save_folder = Path(f"/work/PRTNR/CHUV/DIR/rgottar1/spatial/Cluster_Benchmark/Visualization/010_results")


def get_data_results(dataset_path, data, save_folder):

    data_path = dataset_path / data
    #Get adata object

    ad_list = get_result_anndata(data_path, methods=all_methods)

    for adata in ad_list:
        # adata.obs = adata.obs.T.drop_duplicates().T
        save_dir = save_folder/ data
        save_dir.mkdir(exist_ok=True, parents=True)

        adata.obs = adata.obs.astype(str)
        adata.write_h5ad(save_dir / f"adata_{adata.uns['sample']}.h5ad")
        adata.obs.to_csv(save_dir / f"obs_{adata.uns['sample']}.tsv", sep="\t", index_label="")

# %% run all the dataset results
data_list = ["libd_dlpfc", 
             "visium_breast_cancer_SEDR",
             "visium_chicken_heart", 
             "STARmap_plus",
             "abc_atlas_wmb_thalamus",
             "xenium-mouse-brain-SergioSalas"]

for data in data_list:
    get_data_results(dataset_path=dataset_path, 
                     data=data, 
                     save_folder=save_folder)

# %%
