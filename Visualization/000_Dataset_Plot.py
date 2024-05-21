#%% Load packages
import os
from pathlib import Path
import matplotlib.pyplot as plt
import scanpy as sc

import warnings
warnings.filterwarnings("ignore")

#%% define function to get results
def get_anndata(folder, samples, qc=False):
    # Untested template
    import anndata as ad
    import pandas as pd
    import scipy as sp
    from pathlib import Path

    ad_list = []
    for sample in samples:
        sample_path = folder / sample
        sample_path_qc = Path(sample_path) / "qc" if qc else Path(sample_path)

        X = sp.io.mmread(sample_path_qc / "counts.mtx")
        if sp.sparse.issparse(X):
            X = X.tocsr()
        observations = pd.read_table(sample_path_qc / "observations.tsv", index_col=0)
        labels = pd.read_table(sample_path / "labels.tsv", index_col=0)
        observations = observations.join(labels)
        features = pd.read_table(sample_path_qc / "features.tsv", index_col=0)
        coordinates = (
            pd.read_table(sample_path_qc / "coordinates.tsv", index_col=0)
            .loc[observations.index, :]
            .to_numpy()
        )

        adata = ad.AnnData(
            X=X, obs=observations, var=features, 
            obsm={"spatial": coordinates}, uns={"sample":sample}
        )
        ad_list.append(adata)

    return ad_list

#%% Define methods and datasets
dataset_path = Path("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Cluster_Benchmark/data/")

# data = "libd_dlpfc"
# data = "visium_breast_cancer_SEDR"
# data = "visium_chicken_heart"
# data = "STARmap_plus"
# data = "abc_atlas_wmb_thalamus"
# data = "xenium-mouse-brain-SergioSalas"
data = "stereoseq_developing_Drosophila_embryos_larvae"

data_path = dataset_path / data

samples = [m for m in os.listdir(data_path) if os.path.isdir(os.path.join(data_path, m)) and not m.startswith(".")]

#%% Get adata object
ad_list = get_anndata(data_path, samples, qc=True)

#%% plot ground_truth label used in each sample
def plot_label_spatial(ad_list, spot_size, data, nrow=None, ncol=None, save=True, show=True):
    import math

    nrow = math.ceil(math.sqrt(len(ad_list))) if nrow is None else nrow
    ncol = math.ceil(len(ad_list) / nrow) if ncol is None else ncol
    fig, axs = plt.subplots(nrows=nrow, ncols=ncol, figsize=(25, 13))
    axf = axs.flatten()

    save_dir = Path(f"/work/PRTNR/CHUV/DIR/rgottar1/spatial/Cluster_Benchmark/Visualization/000_Label_Plots")
    save_dir.mkdir(exist_ok=True, parents=True)

    for i, adata in enumerate(ad_list):
        # adata.obs = adata.obs.T.drop_duplicates().T
        sc.pl.spatial(adata, 
                      color="label", 
                      spot_size = spot_size, 
                      title = adata.uns["sample"],
                      show=False,
                      ax = axf[i])
    fig.suptitle(f"Ground Truth labels for {data}")
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    if save:
        plt.savefig(save_dir / f"{data}_label.png", bbox_inches="tight")

    if show:
        plt.show()
    else:
        return fig

plot_label_spatial(ad_list, spot_size=0.5, data=data)

# %%
