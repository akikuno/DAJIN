import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import logomaker as lm

raw_seqs1 = pd.read_csv('.tmp_/tmp1_Fw.fa', sep="\t", header=None)
raw_seqs2 = pd.read_csv('.tmp_/tmp2_Fw.fa', sep="\t", header=None)
expected = pd.read_csv('.tmp_/mutation_Fw.fa', sep="\t", header=None)

seq1 = raw_seqs1[~raw_seqs1[0].str.startswith(">")].squeeze()
seq2 = raw_seqs2[~raw_seqs2[0].str.startswith(">")].squeeze()
expected2 = expected[~expected[0].str.startswith(">")].squeeze()
type(pd.Series(expected2))

counts_mat1 = lm.alignment_to_matrix(seq1, to_type="probability")
counts_mat2 = lm.alignment_to_matrix(seq2, to_type="probability")
counts_expected = lm.alignment_to_matrix(pd.Series(expected2))
counts_mat1.head()

plt.rcParams["font.size"] = 20
plt.rcParams["font.family"] = "Arial"
fig = plt.figure(figsize=(18, 7))
ax1 = fig.add_subplot(3, 1, 1)
ax2 = fig.add_subplot(3, 1, 2)
ax3 = fig.add_subplot(3, 1, 3)
fig.suptitle("Barcode04", fontsize=30)
ax3.set_xlabel("Position", fontsize=25)
ax2.set_ylabel("Probability", fontsize=25)
#
ax1.set_title("Expected joint sequence")
ax2.set_title("Probable intact sequence")
ax3.set_title("Probable non-intact sequence")
plt.tight_layout()
plt.subplots_adjust(top=0.85)


logo1 = lm.Logo(counts_expected, font_name='monospace', ax=ax1,
                color_scheme="colorblind_safe", width=0.9, vpad=0.1)
logo2 = lm.Logo(counts_mat1, font_name='monospace', ax=ax2,
                color_scheme="colorblind_safe", width=0.9, vpad=0.1)
logo3 = lm.Logo(counts_mat2, font_name='monospace', ax=ax3,
                color_scheme="colorblind_safe", width=0.9, vpad=0.1)
logo1.style_xticks(spacing=5)
logo2.style_xticks(spacing=5)
logo3.style_xticks(spacing=5)

plt.savefig("figure.png")
