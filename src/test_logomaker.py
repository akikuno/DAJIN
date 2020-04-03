import sys
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
plt.style.use('ggplot')
plt.rcParams.update({'font.size': 15})
plt.tight_layout()
import seaborn as sns
sns.set(style='ticks', context='talk')

import logomaker as lm

# TEST <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
target_1 = "CCCACAGCTAGGGCCGCATAACTTCGTATAATGTATGCTATACGAAGTTATAGCTTGATATCGAATTGGCGCGCAGTGGCTG"
fasta_intact = pd.read_csv(f".DAJIN_temp/seqlogo/tmp_lalign_barcode14_target",
                        sep=" ", header=None, names=["seq"])
# seq_intact = fasta_intact[~fasta_intact.fa.str.startswith(">")].squeeze()
#seq_max = fasta_intact["seq"].str.len().max()
#fasta_target = target_1.ljust(seq_max, "-")
# fasta_intact["seq"] = fasta_intact["seq"].str.ljust(seq_max, "-")

counts_expected = lm.alignment_to_matrix(pd.Series(target_1),
                                       to_type="counts",
                                       characters_to_ignore='.-X')
counts_intact = lm.alignment_to_matrix(fasta_intact["seq"],
                                       to_type="counts",
                                       characters_to_ignore='.-X')
#plt.subplots_adjust(top=0.85, hspace=1, wspace=1)
fig, ax = plt.subplots(2,1,figsize=(40, 5))

logo_exp = lm.Logo(counts_expected, font_name='monospace', ax=ax[0],
                color_scheme="colorblind_safe", width=0.9, vpad=0.1)
logo_exp.style_xticks(spacing=5)
logo1 = lm.Logo(counts_intact, font_name='monospace', ax=ax[1],
                color_scheme="colorblind_safe", width=0.9, vpad=0.1)
logo1.style_xticks(spacing=5)
plt.savefig(f"test_seqlogo.png")

#
# TEST <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<






# ============================================================================
# # arguments
# ============================================================================
args = sys.argv
fasta_expected = pd.read_csv(args[1], sep="\t", header=None)
fasta_intact = pd.read_csv(args[2], sep="\t", header=None)
fasta_nonintact = pd.read_csv(args[3], sep="\t", header=None)
alignment = pd.read_csv(args[4], sep="\t", header=None)
intact_ratio = pd.read_csv(args[5], sep="\t", header=None)
#
alignment = alignment.T
alignment.columns = ["all", "aligned"]
alignment["non-aligned"] = alignment["all"] - alignment["aligned"]
alignment["per_align"] = alignment["aligned"]/alignment["all"]*100
alignment["per_non-align"] = alignment["non-aligned"]/alignment["all"]*100
#
intact_ratio = intact_ratio.T
intact_ratio.columns = ["intact", "non-intact"]
intact_all = intact_ratio["intact"] + intact_ratio["non-intact"]
intact_ratio["per_intact"] = intact_ratio["intact"]/intact_all*100
intact_ratio["per_non-intact"] = intact_ratio["non-intact"]/intact_all*100
# ============================================================================
# # output figure name
# ============================================================================
output_figname = re.sub(r".*_intact_", "", args[2])
output_figname = re.sub(".fa", "", output_figname)

# TEST =====================================
# barcode = "barcode02"
# fasta_expected = pd.read_csv(".tmp_/mutation.fa", sep="\t", header=None)
# fasta_intact = pd.read_csv(
#     f".tmp_/lalign_intact_{barcode}.fa", sep="\t", header=None)
# fasta_nonintact = pd.read_csv(
#     f".tmp_/lalign_nonintact_{barcode}.fa", sep="\t", header=None)
# alignment = pd.read_csv(f".tmp_/numseq_alignment_{barcode}",
#                         sep="\t", header=None, names=[""])
# intact_ratio = pd.read_csv(f".tmp_/numseq_intact_{barcode}",
#                            sep="\t", header=None, names=[""])
# #
# alignment = alignment.T
# alignment.columns = ["all", "aligned"]
# alignment["non-aligned"] = alignment["all"] - alignment["aligned"]
# alignment["per_align"] = alignment["aligned"]/alignment["all"]*100
# alignment["per_non-align"] = alignment["non-aligned"]/alignment["all"]*100

# intact_ratio = intact_ratio.T
# intact_ratio.columns = ["intact", "non-intact"]
# intact_all = intact_ratio["intact"] + intact_ratio["non-intact"]
# intact_ratio["per_intact"] = intact_ratio["intact"]/intact_all*100
# intact_ratio["per_non-intact"] = intact_ratio["non-intact"]/intact_all*100

# output_figname = re.sub(r".*_intact_", "", f".tmp_/lalign_intact_{barcode}.fa")
# output_figname = re.sub(".fa", "", output_figname)
# TEST <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# # TEST =====================================
# barcode = "barcode02"
# site = "right"
# fasta_expected = pd.read_csv(
#     f".tmp_/cutting_sites_{site}.fa", sep="\t", header=None)
# fasta_intact = pd.read_csv(
#     f".tmp_/lalign_intact_{barcode}_{site}.fa", sep="\t", header=None)
# fasta_nonintact = pd.read_csv(
#     f".tmp_/lalign_nonintact_{barcode}_{site}.fa", sep="\t", header=None)
# alignment = pd.read_csv(f".tmp_/numseq_alignment_{barcode}_{site}",
#                         sep="\t", header=None, names=[""])
# intact_ratio = pd.read_csv(f".tmp_/numseq_intact_{barcode}_{site}",
#                            sep="\t", header=None, names=[""])
# #
# alignment = alignment.T
# alignment.columns = ["all", "aligned"]
# alignment["non-aligned"] = alignment["all"] - alignment["aligned"]
# alignment["per_align"] = alignment["aligned"]/alignment["all"]*100
# alignment["per_non-align"] = alignment["non-aligned"]/alignment["all"]*100

# intact_ratio = intact_ratio.T
# intact_ratio.columns = ["intact", "non-intact"]
# intact_all = intact_ratio["intact"] + intact_ratio["non-intact"]
# intact_ratio["per_intact"] = intact_ratio["intact"]/intact_all*100
# intact_ratio["per_non-intact"] = intact_ratio["non-intact"]/intact_all*100

# output_figname = re.sub(
#     r".*_intact_", "", f".tmp_/lalign_intact_{barcode}_{site}.fa")
# output_figname = re.sub(".fa", "", output_figname)
# title_figname = re.sub(
#     r"_", " ", output_figname)

# ============================================================================
# # create position weight matrix
# ============================================================================

seq_expected = fasta_expected[~fasta_expected[0].str.startswith(">")].squeeze()
seq_intact = fasta_intact[~fasta_intact.fa.str.startswith(">")].squeeze()
seq_nonintact = fasta_nonintact[~fasta_nonintact[0].str.startswith(
    ">")].squeeze()

counts_seq_expected = lm.alignment_to_matrix(
    pd.Series(seq_expected), characters_to_ignore='.-X')
counts_intact = lm.alignment_to_matrix(pd.Series(seq_intact),
                                       to_type="counts",
                                       characters_to_ignore='.-X')
counts_nonintact = lm.alignment_to_matrix(pd.Series(seq_nonintact),
                                          to_type="counts",
                                          characters_to_ignore='.-X')

# ============================================================================
# # matplotlib #? ----------------------------------------------
# ============================================================================

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Params
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
plt.rcParams["font.size"] = 18
plt.rcParams["font.family"] = "Arial"
plt.tight_layout()
#plt.subplots_adjust(top=0.85, hspace=1, wspace=1)
fig = plt.figure(figsize=(18, 10))
gs_master = GridSpec(nrows=3, ncols=3, width_ratios=[1, 1, 8],
                     hspace=0.1, wspace=0.4)

gs_1 = GridSpecFromSubplotSpec(
    nrows=3, ncols=1, subplot_spec=gs_master[0:2, 0], wspace=0.4)
ax1 = fig.add_subplot(gs_1[:, :])

gs_2 = GridSpecFromSubplotSpec(
    nrows=3, ncols=1, subplot_spec=gs_master[0:2, 1], wspace=0.4)
ax2 = fig.add_subplot(gs_2[:, :])

gs_3_4_5 = GridSpecFromSubplotSpec(
    nrows=3, ncols=1, subplot_spec=gs_master[0:2, 2], hspace=1, wspace=1)
ax3 = fig.add_subplot(gs_3_4_5[0, :])
ax4 = fig.add_subplot(gs_3_4_5[1, :])
ax5 = fig.add_subplot(gs_3_4_5[2, :])

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Title, labels
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

fig.suptitle(title_figname, fontsize=35)
ax1.set_ylabel("The number of reads", fontsize=20)
ax2.set_ylabel("The number of reads", fontsize=20)
ax3.set_title("Expected sequence")
ax4.set_title("Probable intact sequence")
ax4.set_ylabel("Counts", fontsize=20)
ax5.set_title("Probable non-intact sequence")
ax5.set_xlabel("Position", fontsize=20)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Create figures
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Stacked plot
alignment[["aligned", "non-aligned"]].plot.bar(
    stacked=False, ax=ax1,
    color=["darkorange", "steelblue"]).legend(
    bbox_to_anchor=(0.5, -0.2), loc="lower center", framealpha=1)

intact_ratio[["intact", "non-intact"]].plot.bar(
    stacked=False, ax=ax2,
    color=["darkorange", "steelblue"]).legend(
    bbox_to_anchor=(0.5, -0.2), loc="lower center", framealpha=1)

# SEQ LOGO
logo1 = lm.Logo(counts_seq_expected, font_name='monospace', ax=ax3,
                color_scheme="colorblind_safe", width=0.9, vpad=0.1)
logo1.style_xticks(spacing=5)
#
if int(intact_ratio["intact"]) > int(alignment["all"] * 0.05):
    logo2 = lm.Logo(counts_intact, font_name='monospace', ax=ax4,
                    color_scheme="colorblind_safe", width=0.9, vpad=0.1)
    logo2.style_xticks(spacing=5)
else:
    ax4.text(0.5, 0.5,
             "Too few reads to display",
             horizontalalignment='center',
             verticalalignment='center',
             size=30, color="black")

#
if int(intact_ratio["non-intact"]) > int(alignment["all"] * 0.05):
    logo3 = lm.Logo(counts_nonintact, font_name='monospace', ax=ax5,
                    color_scheme="colorblind_safe", width=0.9, vpad=0.1)
    logo3.style_xticks(spacing=5)
else:
    ax5.text(0.5, 0.5,
             "Too few reads to display",
             horizontalalignment='center',
             verticalalignment='center',
             size=30, color="black")


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Save figure
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
plt.savefig(f"test_seqlogo_{output_figname}.png")
