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
df_predicted = pd.read_csv(".tmp_/prediction_result.txt",
                           sep="\t")
barcode_list = df_mutation.barcodeID.unique()

# ? ==========================================================

df_read_size = df_predicted.groupby("barcodeID").size()
df_mutation_size = df_mutation.groupby("barcodeID").size()

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax1.grid(True)
ax2.grid(False)
# ax1.set_axisbelow(True)

color_1 = plt.cm.Set1.colors[1]
color_1 = "#696969"
color_2 = plt.cm.Set1.colors[0]

df_read_size.plot(ax=ax1, kind='bar', stacked=False, rot=0,
                  color=color_1, label="total read numbers")
df_mutation_size.plot(ax=ax2, kind='bar', stacked=False,
                      rot=0, color=color_2, label="total flox numbers")

ax1.yaxis.set_major_locator(MaxNLocator(nbins=5))
ax2.yaxis.set_major_locator(MaxNLocator(nbins=5))

# ax2.spines['left'].set_color(color_1)
# ax2.spines['right'].set_color(color_2)

# ax1.tick_params(axis='y', colors=color_1)
#ax2.tick_params(axis='y', colors=color_2)

handler1, label1 = ax1.get_legend_handles_labels()
handler2, label2 = ax2.get_legend_handles_labels()

ax1.legend(handler1 + handler2, label1 + label2, loc=2, borderaxespad=0.)

read_max = 1.2 * max(df_read_size)
ax1.set_ylim([0, read_max])
ax2.set_ylim([0,  read_max])
ax1.set_xticklabels(ax1.get_xticklabels(), rotation=90)
ax2.set_xticklabels(ax2.get_xticklabels(), rotation=90)

fig_name = "fullmatchsearch_read_numbers"
for fig_type in fig_types:
    plt.savefig(f"{output_dir}/{fig_type}/{fig_name}.{fig_type}",
                dpi=350, bbox_inches="tight")
# -


#! The percentage of flox fullmatch reads ----------------------------

df_stacked = np.zeros(
    [df_mutation.barcodeID.value_counts().max(), len(barcode_list)])
df_stacked = pd.DataFrame(df_stacked, columns=barcode_list)

for i in barcode_list:
    df_stacked[i] = df_mutation[df_mutation.barcodeID ==
                                i]["mutation"].reset_index(drop=True)

# +

counts = df_stacked.apply(lambda x: x.dropna(
).value_counts() / len(x.dropna())).transpose()

counts_flox = counts["matchedx2"] * 100

sns.set_style("ticks", {"font": "Arial"})
plt.style.use('seaborn-ticks')
sns.set_palette("Paired")
fig = plt.figure()
ax = fig.add_subplot(111)

counts_flox.plot(ax=ax, kind='bar', stacked=True, rot=0, color="#FF4500")
ax.set_ylabel("Percentage of flox fullmatch reads")
ax.yaxis.grid(True)
ax.set_axisbelow(True)
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
#ax.legend(bbox_to_anchor=(1, 1, 0.1, 0))

# figure ----------------------------
fig_name = "mutation_percentage"
for fig_type in fig_types:
    plt.savefig(f"{output_dir}/{fig_type}/{fig_name}.{fig_type}",
                dpi=350, bbox_inches="tight")
# -


#! ## The percentage of flox fullmatched reads --------------------

df_stacked = np.zeros(
    [df_mutation.barcodeID.value_counts().max(), len(barcode_list)])
df_stacked = pd.DataFrame(df_stacked, columns=barcode_list)

for i in barcode_list:
    df_stacked[i] = df_mutation[df_mutation.barcodeID ==
                                i]["mutation"].reset_index(drop=True)

# +

counts = df_stacked.apply(lambda x: x.dropna(
).value_counts() / len(x.dropna())).transpose()

sns.set_style("ticks", {"font": "Arial"})
plt.style.use('seaborn-ticks')
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111)

# colorlist = ["#FFB6C1", "#FF4500", "#DCDCDC"]
colorlist = sns.color_palette("colorblind", len(counts.columns)-1)
colorlist.append("#DDDDDD")
colorlist
# counts.plot(ax=ax, kind='bar', stacked=True, rot=0, color=colorlist)
# vals = ax.get_yticks()
# ax.set_yticklabels(['{:3.2f}%'.format(x*100) for x in vals])
# ax.yaxis.grid(True)
# ax.set_axisbelow(True)
# ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
# ax.legend(bbox_to_anchor=(1, 1, 0.1, 0))

counts.iloc[::-1].plot(ax=ax, kind='barh', stacked=True, rot=0,
                       color=colorlist,
                       edgecolor='#000000', width=1)  # "#f87f73"

ax.set_xlabel("Percentage of predicted allele type")
ax.set_xticklabels(['{:3.0f}%'.format(x*100) for x in ax.get_xticks()])
ax.yaxis.grid(True)
ax.set_axisbelow(True)
ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
ax.legend(bbox_to_anchor=(1, 1, 0.1, 0))

# figure ----------------------------
fig_name = "mutation_fullmatch"
for fig_type in ["svg", "png"]:
    plt.savefig(f"{output_dir}/{fig_type}/{fig_name}.{fig_type}",
                dpi=350, bbox_inches="tight")

# +

#! ## Numbers of predicted flox reads and all reads ----------------

df_mutation.head()

df_mutation_size = df_mutation.groupby("barcodeID").size()

# +
sns.set_style("ticks", {"font": "Arial"})
plt.style.use('seaborn-ticks')
sns.set_palette("Paired")

fig = plt.figure()
ax = fig.add_subplot(111)

df_mutation_size.plot(ax=ax, kind='bar', stacked=True, rot=0, color="#333631")
# vals = ax.get_yticks()
# ax.set_yticklabels(['{:3.2f}%'.format(x*100) for x in vals])
ax.set_ylabel("Number of predicted flox reads")
ax.yaxis.grid(True)
ax.set_axisbelow(True)
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
#ax.legend(bbox_to_anchor=(1, 1, 0.1, 0))

# figure ----------------------------
fig_name = "total_flox_numbers"
for fig_type in fig_types:
    plt.savefig(f"{output_dir}/{fig_type}/{fig_name}.{fig_type}",
                dpi=350, bbox_inches="tight")
# -
# -

#! The number of total analyzed reads ------------------------------------

df_read_size = df_predicted.groupby("barcodeID").size()

# +
sns.set_style("ticks", {"font": "Arial"})
plt.style.use('seaborn-ticks')
sns.set_palette("Paired")

fig = plt.figure()
ax = fig.add_subplot(111)

df_read_size.plot(ax=ax, kind='bar', stacked=True, rot=0, color="#333631")
# vals = ax.get_yticks()
# ax.set_yticklabels(['{:3.2f}%'.format(x*100) for x in vals])
ax.set_ylabel("Number of total analyzed reads")
ax.yaxis.grid(True)
ax.set_axisbelow(True)
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
#ax.legend(bbox_to_anchor=(1, 1, 0.1, 0))

# figure ----------------------------
fig_name = "read_numbers"
for fig_type in fig_types:
    plt.savefig(f"{output_dir}/{fig_type}/{fig_name}.{fig_type}",
                dpi=350, bbox_inches="tight")
# -
