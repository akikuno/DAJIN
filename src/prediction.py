from keras.models import Model
from keras import Model
from keras.models import load_model
import keras
from matplotlib.axes._axes import _log as matplotlib_axes_logger
import sys
import itertools
from sklearn.metrics import confusion_matrix
from keras.models import Sequential
from keras import regularizers
from keras.layers import Conv1D, Dense, MaxPooling1D, Flatten, Dropout, Activation
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from keras.utils import np_utils
from sklearn.model_selection import train_test_split
from functools import partial
import os
import numpy as np
import pandas as pd
import re
from sklearn.preprocessing import OneHotEncoder
from sklearn.preprocessing import LabelEncoder
from numpy import argmax
from numpy import array
from keras.backend import tensorflow_backend
from tensorflow.python.client import device_lib
import tensorflow as tf
import warnings
from tqdm import tqdm
warnings.filterwarnings('ignore')

# ====================================
# # Import data
# ====================================

args = sys.argv
file_name = args[1]
# file_name = "data_for_ml/sequence_MIDS.txt.gz"

# df_anomaly = pd.read_csv(args[2], header=None, sep='\t')
df_anomaly = pd.read_csv(".tmp_/anomaly_classification_revised.txt",
                         header=None, sep='\t')
labels_index = pd.read_csv(
    ".tmp_/anomaly_classification_labels.txt",
    header=None, sep='\t')
# labels_id = df_fasta[df_fasta.iloc[:, 0].str.startswith(
#     ">")].squeeze().str.strip(">").reset_index(drop=True)

fig_dirs = ["results/figures/png", "results/figures/svg"]

output_npz = file_name.replace(".txt.gz", ".npz").replace(
    "data_for_ml/", "data_for_ml/model/")
output_figure = file_name.replace(".txt.gz", "").replace("data_for_ml/", "")
output_model = file_name.replace(".txt.gz", "").replace(
    "data_for_ml", "data_for_ml/model")

# ====================================
# # Load One-hot matrix
# ====================================
np.load = partial(np.load, allow_pickle=True)
npz = np.load(output_npz)
X_real = npz["X_real"]

# ====================================
# # Load model
# ====================================

tf.get_logger().setLevel('INFO')
model = load_model(output_model + '.h5')

# ====================================
# # Prediction
# ====================================

print("Predict labels...")
predict = model.predict(X_real, verbose=1, batch_size=1)
# predict = model.predict(X_real[0:1000], verbose=1, batch_size=1)

predict = np.argmax(predict, axis=1)

df_predict = pd.Series(predict, dtype="str") + "_"

for i, j in enumerate(labels_index.squeeze().str.replace("_simulated", "")):
    df_predict = df_predict.str.replace(str(i)+"_", j)

df_result = pd.DataFrame({"barcodeID": df_anomaly.iloc[:, 0],
                          "seqID": df_anomaly.iloc[:, 1],
                          "predict": df_predict,
                          "anomaly": df_anomaly.iloc[:, 2]})

df_result.predict = df_result.predict.mask(
    df_result.anomaly.str.contains("Abnormal"), df_result.anomaly)

del df_result["anomaly"]
# df_result = df_result.head(1000)

# ====================================
# ## Output result
# ====================================

df_result.to_csv('.tmp_/prediction_result.txt', sep='\t', index=False)

# ====================================
# ## Visualization of allele profile
# ====================================

barcode_list = df_result.barcodeID.unique()
df_stacked = np.zeros(
    [df_result.barcodeID.value_counts().max(), len(barcode_list)])
df_stacked = pd.DataFrame(df_stacked, columns=barcode_list)
for i in barcode_list:
    df_stacked[i] = df_result[df_result.barcodeID
                              == i]["predict"].reset_index(drop=True)

# ## Plot figures
colorlist = ["#FF4500", "#D3D3D3", "#ADD8E6"]  # #88CCEE #0072B2
colorlist = ["#FF4500", "#D3D3D3", "#FFF9B0", "#ADD8E6"]
colorlist.extend(list(sns.color_palette("Accent", 24).as_hex()))

counts = df_stacked.apply(lambda x: x.dropna(
).value_counts() / len(x.dropna())).transpose()

tmp0 = pd.Series(["target", "wt"])
tmp1 = counts.columns[counts.columns.str.startswith(("Abnormal"))].values
tmp10 = pd.concat([tmp0, pd.Series(tmp1)])

tmp2 = counts.drop(tmp10, axis=1).columns.values
tmp012 = pd.concat([tmp10, pd.Series(tmp2)])

counts = counts[tmp012]

# counts = counts[np.concatenate([tmp0, tmp1, tmp2])]
sns.set_style("ticks", {"font": "Arial"})
plt.style.use('seaborn-pastel')
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111)
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
fig_name = "prediction_result"
for fig_dir in fig_dirs:
    fig_type = re.sub(".*/", "", fig_dir)
    plt.savefig(f"{fig_dir}/{output_figure}_{fig_name}.{fig_type}",
                dpi=350, bbox_inches="tight")
