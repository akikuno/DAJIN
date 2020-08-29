import os

os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"
os.environ["TF_FORCE_GPU_ALLOW_GROWTH"] = "true"

import sys
import random as rn
import numpy as np
import pandas as pd

import pickle

from sklearn.preprocessing import LabelEncoder, LabelBinarizer

import tensorflow as tf
from tensorflow.keras.models import Model

# set random seed
os.environ['PYTHONHASHSEED'] = '0'
np.random.seed(42)
rn.seed(12345)
tf.random.set_seed(1234)

################################################################################
#! I/O naming
################################################################################

# ===========================================================
# ? TEST auguments
# ===========================================================

# file_name = ".DAJIN_temp/data/MIDS_barcode11_wt"
# mutation_type = "S"
# threads = 12

# ===========================================================
# ? Auguments
# ===========================================================

args = sys.argv
file_name = args[1]
mutation_type = args[2]
threads = int(args[3])

if mutation_type == "":
    raise ValueError("mutation_type is empty")

if threads == "":
    import multiprocessing
    threads = multiprocessing.cpu_count() // 2


# ===========================================================
# ? Input
# ===========================================================

df = pd.read_csv(file_name, header=None, sep="\t")
df.columns = ["seqID", "seq", "barcodeID"]
df.seq = "MIDS=" + df.seq

# ===========================================================
# ? Encording Function
# ===========================================================


def label_encorde_seq(seq):
    label_seq = (
        seq.apply(list).apply(LabelEncoder().fit_transform).apply(pd.Series).to_numpy()
    )
    return label_seq


def onehot_encode_seq(seq):
    onehot_seq = np.apply_along_axis(
        LabelBinarizer().fit_transform, 1, label_encorde_seq(seq)
    ).astype(np.uint8)
    return onehot_seq


X_real = onehot_encode_seq(df["seq"])

del df["seq"]  # <<<

################################################################################
#! load trained models
################################################################################

labels_index = np.load(".DAJIN_temp/data/labels_index.npy", allow_pickle=True)
model = tf.keras.models.load_model(".DAJIN_temp/data/model.h5")
clf = pickle.load(open(".DAJIN_temp/data/model_lof.sav", "rb"))

################################################################################
#! Novelity (Anomaly) detection
################################################################################
# ===========================================================
# ? Extract layer
# ===========================================================

model_ = Model(model.get_layer(index=0).input, model.get_layer(index=-2).output)
predict_vector = model_.predict(X_real, verbose=0, batch_size=32)

# ===========================================================
# ? LocalOutlierFactor
# ===========================================================

outliers = clf.predict(predict_vector)
outliers = np.where(outliers == 1, "normal", "abnormal")

df["outliers"] = outliers

df.groupby("barcodeID").outliers.value_counts() #!<<<<<


################################################################################
#! Prediction
################################################################################

predict = model.predict(X_real.astype("float16"), verbose=0, batch_size=32)
predict = np.argmax(predict, axis=1)

# del X_real  # <<<

df["prediction"] = predict
df["prediction"].mask(df["outliers"] == "abnormal", "abnormal", inplace=True)

del df["outliers"]  # <<<

for index, label in enumerate(labels_index):
    label = label.replace("_simulated", "")
    df["prediction"].mask(df["prediction"] == index, label, inplace=True)

# ---------------------------------------
# * In the case of a point mutation,
# * the reads named "wt_ins" and "wt_del"
# * should be treated as "abnormal".
# ---------------------------------------
if mutation_type == "S":
    df["prediction"].mask(
        df["prediction"].str.contains("wt_"), "abnormal", inplace=True
    )

df.groupby("barcodeID").prediction.value_counts() #!<<<<<

################################################################################
#! Save file
################################################################################

df.to_csv(
    ".DAJIN_temp/data/tmp_DAJIN_MIDS_prediction_result.txt",
    header=False,
    index=False,
    sep="\t",
    mode="a",
)

