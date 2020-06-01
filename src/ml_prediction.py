import os
os.environ['TF_CPP_MIN_LOG_LEVEL']='3'
import sys
import re
from functools import partial
import warnings
warnings.filterwarnings("ignore")
warnings.simplefilter(action='ignore', category=FutureWarning)

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm

from tensorflow.keras.models import load_model


# ====================================
# # Import data
# ====================================

args = sys.argv
file_name = args[1]
# file_name = ".DAJIN_temp/data/DAJIN_MIDS.txt"

df_anomaly = pd.read_csv(".DAJIN_temp/data/DAJIN_anomaly_classification.txt",
                         header=None, sep='\t')
labels_index = pd.read_csv(
    ".DAJIN_temp/data/DAJIN_anomaly_classification_labels.txt",
    header=None, sep='\t')
labels_index.columns = ["label"]

output_npz = file_name.replace(".txt", ".npz")
output_model = file_name.replace(".txt", ".h5")

# fig_dirs = ["results/figures/png", "results/figures/svg"]
# output_figure = file_name.replace(".txt", "")

# ====================================
# # Load One-hot matrix
# ====================================
# np.load = partial(np.load, allow_pickle=True)
npz = np.load(output_npz)
X_real = npz["X_real"]

# ====================================
# # Load model
# ====================================

model = load_model(output_model)

# ====================================
# # Prediction
# ====================================
print("Predict allele types...")
iter_ = 1000
predict = np.zeros(X_real.shape[0], dtype="uint8")
for i in tqdm(range(0, X_real.shape[0], iter_)):
    predict_ = model.predict(X_real[i: i + iter_, :, :].astype("float16"),
                             verbose=0, batch_size=32)
    predict[i:i+iter_] = np.argmax(predict_, axis=1)

df_predict = pd.Series(predict, dtype="str") + "_"

for i, j in enumerate(labels_index["label"].str.replace("_simulated.*$", "")):
    df_predict = df_predict.str.replace(str(i)+"_", j)

df_result = pd.DataFrame({"barcodeID": df_anomaly.iloc[:, 0],
                          "seqID": df_anomaly.iloc[:, 1],
                          "predict": df_predict,
                          "anomaly": df_anomaly.iloc[:, 2]})

df_result.predict = df_result.predict.mask(
    df_result.anomaly.str.contains("abnormal"), df_result.anomaly)

del df_result["anomaly"]
# df_result = df_result.head(1000)

# ====================================
# ## Output result
# ====================================

df_result.to_csv('.DAJIN_temp/data/DAJIN_MIDS_prediction_result.txt',
                 sep='\t', index=False, header=False)
                 
# ====================================
# ## Visualization of allele profile
# ====================================
# data = df_result.sort_values(by="barcodeID")
# barcode_list = data.barcodeID.unique()
# df_stacked = np.zeros(
#     [data.barcodeID.value_counts().max(), len(barcode_list)])
# df_stacked = pd.DataFrame(df_stacked, columns=barcode_list)
# for i in barcode_list:
#     df_stacked[i] = data[data.barcodeID == i]["predict"].reset_index(drop=True)

# # ## Plot figures
# colorlist = ["#FF4500", "#DDDDDD", "#ADD8E6"]  # #88CCEE #0072B2 # ccffff
# colorlist = ["#FF4500", "#DDDDDD", "#B0E0E6"]  # #88CCEE #0072B2 #
# # colorlist = ["#FF4500", "#D3D3D3", "#FFF9B0", "#ADD8E6"]
# colorlist.extend(list(sns.color_palette("Accent", 24).as_hex()))

# counts = df_stacked.apply(lambda x: x.dropna(
# ).value_counts() / len(x.dropna())).transpose()

# tmp0 = pd.Series(["target", "wt"])
# tmp1 = counts.columns[counts.columns.str.startswith(("abnormal"))].values
# tmp10 = pd.concat([tmp0, pd.Series(tmp1)])

# tmp2 = counts.drop(tmp10, axis=1).columns.values
# tmp012 = pd.concat([tmp10, pd.Series(tmp2)])

# counts = counts[tmp012]

# # counts = counts[np.concatenate([tmp0, tmp1, tmp2])]
# sns.set_style("ticks", {"font": "Arial"})
# plt.style.use('seaborn-pastel')
# fig = plt.figure(figsize=(10, 10))
# ax = fig.add_subplot(111)
# counts.iloc[::-1].plot(ax=ax, kind='barh', stacked=True, rot=0,
#                        color=colorlist,
#                        edgecolor='#000000', width=1)  # "#f87f73"

# ax.set_xlabel("Percentage of predicted allele type")
# ax.set_xticklabels(['{:3.0f}%'.format(x*100) for x in ax.get_xticks()])
# ax.yaxis.grid(True)
# ax.set_axisbelow(True)
# ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
# ax.legend(bbox_to_anchor=(1, 1, 0.1, 0))

# # figure ----------------------------
# fig_name = "prediction_result"
# for fig_dir in fig_dirs:
#     fig_type = re.sub(".*/", "", fig_dir)
#     plt.savefig(f"{fig_dir}/{output_figure}_{fig_name}.{fig_type}",
#                 dpi=350, bbox_inches="tight")


