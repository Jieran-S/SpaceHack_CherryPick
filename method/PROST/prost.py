#!/usr/bin/env python

# Author_and_contribution: Niklas Mueller-Boetticher; created script

import argparse

parser = argparse.ArgumentParser(
    description="""SpaGCN (https://www.nature.com/articles/s41592-021-01255-8)
Requirements: Visium or ST. Images can be used."""
)

parser.add_argument(
    "-c", "--coordinates", help="Path to coordinates (as tsv).", required=True
)
parser.add_argument(
    "-m", "--matrix", help="Path to (transformed) counts (as mtx).", required=False
)
parser.add_argument(
    "-f", "--features", help="Path to features (as tsv).", required=True
)
parser.add_argument(
    "-o", "--observations", help="Path to observations (as tsv).", required=True
)
parser.add_argument(
    "-n",
    "--neighbors",
    help="Path to neighbor definitions. Square matrix (not necessarily symmetric) where each row contains the neighbors of this observation (as mtx).",
    required=False,
)
parser.add_argument("-d", "--out_dir", help="Output directory.", required=True)
parser.add_argument(
    "--dim_red",
    help="Reduced dimensionality representation (e.g. PCA).",
    required=False,
)
parser.add_argument("--image", help="Path to H&E staining.", required=False)
parser.add_argument(
    "--n_clusters", help="Number of clusters to return.", required=True, type=int
)
parser.add_argument(
    "--technology",
    help="The technology of the dataset (Visium, ST, ...).",
    required=True,
)
parser.add_argument(
    "--seed", help="Seed to use for random operations.", required=True, type=int
)
parser.add_argument(
    "--config",
    help="Optional config file used to pass additional parameters.",
    required=False,
)

args = parser.parse_args()

from pathlib import Path

out_dir = Path(args.out_dir)

# Output files
label_file = out_dir / "domains.tsv"
embedding_file = out_dir / "embedding.tsv"

n_clusters = args.n_clusters
technology = args.technology
seed = args.seed


## Your code goes here
import json

with open(args.config, "r") as f:
    config = json.load(f)

import PROST
import scipy as sp
import pandas as pd
import random
import scanpy as sc
import torch
import numpy as np

def get_anndata(args):
    import anndata as ad
    import scipy.io
    from PIL import Image

    X = sp.io.mmread(args.matrix)
    if sp.sparse.issparse(X):
        X = X.tocsr()

    observations = pd.read_table(args.observations, index_col=0)
    features = pd.read_table(args.features, index_col=0)

    # Filter
    if "selected" in observations.columns:
        X = X[observations["selected"].to_numpy().nonzero()[0], :]
        observations = observations.loc[lambda df: df["selected"]]
    if "selected" in features.columns:
        X = X[:, features["selected"].to_numpy().nonzero()[0]]
        features = features.loc[lambda df: df["selected"]]

    coordinates = (
        pd.read_table(args.coordinates, index_col=0)
        .loc[observations.index, :]
        .to_numpy()
    )

    adata = ad.AnnData(
        X=X, obs=observations, var=features, obsm={"spatial": coordinates}
    )

    if args.image is not None:
        adata.uns["image"] = np.array(Image.open(args.image))
    else:
        adata.uns["image"] = None

    return adata


adata = get_anndata(args)

# Set seed
random.seed(seed)
torch.manual_seed(seed)
np.random.seed(seed)

platform = "visium" if technology == "Visium" else technology
adata = PROST.prepare_for_PI(adata, platform=platform)
adata = PROST.cal_PI(adata, platform=platform)
PROST.spatial_autocorrelation(adata, k=10, permutations=None)

sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
adata = PROST.feature_selection(adata, by = "prost", n_top_genes = 3000)

PROST.run_PNN(adata,
            platform=platform,
            key_added = "PROST",
            init="mclust",                         
            n_clusters = n_clusters,                        
            lap_filter = 2,                                  
            lr = 0.1,                         
            SEED=seed,                          
            max_epochs = 500,                        
            tol = 5e-3,                        
            post_processing = True,                        
            pp_run_times = 3,
            cuda = False)

## Write output
out_dir.mkdir(parents=True, exist_ok=True)

label_df.columns = ["label"]
label_df.to_csv(label_file, sep="\t", index_label="")
