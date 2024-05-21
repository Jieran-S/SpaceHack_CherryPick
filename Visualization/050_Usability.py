#%% Load packages
import os
import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import marsilea as ma

import warnings
warnings.filterwarnings("ignore")

#%% define function to get results
usability = pd.read_csv("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Cluster_Benchmark/Visualization/050_Usability/Usability_scores.csv", 
                        sep=",", index_col=[0,1])

usability = usability.dropna(axis=1, how="all").apply(pd.to_numeric)

# %% Plot the heatmap of usability
colors = ["#bc3033", "#ff8c00", "#b77961", "#919a77", "#394b61","#fbe9dd", "#391802", "#014431", "#b0b8c3", ]
color_dic = {aspect:color for color, aspect in zip(colors, usability.index.get_level_values(0).unique())}

wes_palette_final = [
    "#710421",  # dark Red
    "#Ff5100",  # International Orange
    "#836247",  # coffee light
    "#69e256",  # Pastel Green
    "#F9C846",  # Bright Mustard Yellow
    "#34568B",  # Dark Blue
    "#8e67d8",  # Medium Purple
    "#1d8af9",  # Dodger Blue
    "#005A51",  # Deep Teal
    "#381756",  # grape Purple
    "#9c9353",  # Corn dark yellow
    "#DC143C",  # Crimson
    "#175620",  # Parsley green
    "#293037",  # Space gray
    "#99D6EA",  # Brighter Sky Blue
    "#C44E55",  # Bright Deep Red
]
Method_color = ma.plotter.Colors(usability.columns, 
                                 palette=wes_palette_final, label="Methods")

# reorder the columns based on column sum
usability = usability[usability.sum().sort_values(ascending=False).index]


h = ma.Heatmap(usability, cmap="PuBu", linewidth=1, width=3.5, height=5.5, linecolor="lightgray")
h.hsplit(labels = usability.index.get_level_values(0), order = usability.index.get_level_values(0).unique())
Aspect_color = ma.plotter.Colors(usability.index.get_level_values(0),
                          palette=color_dic, label="aspects")
h.add_right(Aspect_color, size=.2, pad=.05)
h.add_top(ma.plotter.Numbers(usability.sum(axis=0), show_value=False, color = "#b74923"), pad=.05)
h.add_bottom(ma.plotter.Labels(usability.columns), pad=0.1)
h.add_legends()
h.render()


# %%
