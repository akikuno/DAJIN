#!/usr/bin/env python

# # Environment

# +
import matplotlib.ticker as ticker
from matplotlib.ticker import MaxNLocator
from matplotlib.axes._axes import _log as matplotlib_axes_logger
import sys
import seaborn as sns
import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd
import re
import warnings
warnings.filterwarnings('ignore')


output_dir = "results/figures"
fig_types = ["svg", "png"]

# * Import data from mutation fullmatch search ----------------------------------

df_mutation = pd.read_csv(".tmp_/target_perfectmatch.csv",
                          names=["barcodeID", "mutation"])

barcode_list = df_mutation.barcodeID.unique()
df_barcode_list = pd.DataFrame({"barcodeID": df_mutation.barcodeID.unique()})
df_predicted = pd.read_csv(".tmp_/prediction_result.txt",
                           sep="\t")
df_predicted = pd.merge(df_predicted, df_barcode_list,
                        on="barcodeID", how="inner")

# ======================================
# Stacked barplot of each mutation type
# ======================================

df_stacked = np.zeros(
    [df_mutation.barcodeID.value_counts().max(), len(barcode_list)])
df_stacked = pd.DataFrame(df_stacked, columns=barcode_list)

for i in barcode_list:
    df_stacked[i] = df_mutation[df_mutation.barcodeID ==
                                i]["mutation"].reset_index(drop=True)

counts = df_stacked.apply(lambda x: x.dropna(
).value_counts() / len(x.dropna())).transpose()*100

sns.set_style("ticks", {"font": "Arial"})
plt.style.use('seaborn-ticks')
labels = counts.columns.values
# colorlist = ["#DDDDDD", "red", "darkorange", "gold"]
colorlist = {"No exact match": "#DDDDDD",
             "exact flox": "red",
             "exact left loxP": "darkorange",
             "exact right loxP":"gold"}
colors = []
for i in range(len(labels)):
    colors.append(colorlist[labels[i]])


fig = plt.figure()
ax = fig.add_subplot(111)
counts.plot(ax=ax, kind='bar', stacked=True, rot=0, color=colors)
ax.legend(bbox_to_anchor=(1, 1, 0.1, 0))
ax.set_axisbelow(True)
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
ax.set_ylabel("Percentage of exactly matched reads", fontsize=15)
ax.yaxis.grid(True)
ax.legend(bbox_to_anchor=(1.05, 1), loc = "upper left")

# figure ----------------------------
fig_name = "persentage_of_loxP_intactness"
for fig_type in fig_types:
    plt.savefig(f"{output_dir}/{fig_type}/{fig_name}.{fig_type}",
                dpi=350, bbox_inches="tight")

