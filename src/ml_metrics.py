import os

os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"
os.environ["TF_FORCE_GPU_ALLOW_GROWTH"] = "true"

import sys
import numpy as np
import pandas as pd

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder, LabelBinarizer
from sklearn.neighbors import LocalOutlierFactor

import tensorflow as tf
from tensorflow.keras.layers import Input, Conv1D, Dense, Flatten, MaxPooling1D
from tensorflow.keras.models import Model

################################################################################
#! I/O naming
################################################################################

#==========================================================
#? TEST auguments
#==========================================================

# file_cont = ".DAJIN_temp/data/DAJIN_MIDS_sim.txt"
# file_ab = ".DAJIN_temp/data/DAJIN_MIDS_ab.txt"
# threads = 12
# L2="on"

#==========================================================
#? Auguments
#==========================================================

args = sys.argv
file_cont = args[1]
file_ab = args[2]

if args[3] == "":
    import multiprocessing
    threads = multiprocessing.cpu_count() // 2
else:
    threads = int(args[3])

#==========================================================
#? Input
#==========================================================

df_sim = pd.read_csv(file_cont, header=None, sep="\t")
df_sim.columns = ["seqID", "seq", "barcodeID"]

if "MIDS" in file_cont:
    df_sim.seq = "MIDS=" + df_sim.seq
else:
    df_sim.seq = "ACGT=" + df_sim.seq

#==========================================================
#? Output
#==========================================================



################################################################################
#! Training model
################################################################################

#==========================================================
#? Encording Function
#==========================================================


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


#==========================================================
#? Train test split
#==========================================================

labels, labels_index = pd.factorize(df_sim.barcodeID)
labels = tf.keras.utils.to_categorical(labels)

X_train, X_val, Y_train, Y_val = train_test_split(
    onehot_encode_seq(df_sim.seq), labels, test_size=0.2, shuffle=True
)

#===========================================================
#? Model construction
#===========================================================

inputs = Input(shape = (X_train.shape[1], X_train.shape[2]))
init_kernel_size = int(256)

x = Conv1D(
        filters=32,
        kernel_size=init_kernel_size,
        activation="relu",
        name="1st_Conv1D",
        padding="same",
    )(inputs)
x = MaxPooling1D(pool_size=12, padding="same", name="1st_MaxPooling1D")(x)

x = Conv1D(
        filters=32,
        kernel_size=int(init_kernel_size / 2),
        activation="relu",
        name="2nd_Conv1D",
        padding="same",
    )(x)
x = MaxPooling1D(pool_size=6, padding="same", name="2nd_MaxPooling1D")(x)

x = Conv1D(
        filters=32,
        kernel_size=int(init_kernel_size / 4),
        activation="relu",
        name="3rd_Conv1D",
        padding="same",
    )(x)
x = MaxPooling1D(pool_size=3, padding="same", name="3rd_MaxPooling1D")(x)

x = Flatten(name="flatten")(x)

x = Dense(32, activation="relu", name="1st_FC")(x)

predictions = Dense(len(labels_index), activation="softmax", name="softmax")(x)

model = Model(inputs = inputs, outputs = predictions)

model.compile(optimizer="adam", loss="categorical_crossentropy", metrics=["accuracy"])

# model.summary()
#===========================================================
#? Training
#===========================================================

model.fit(
    X_train,
    Y_train,
    epochs=20,
    verbose=1,
    batch_size=32,
    validation_data=(X_val, Y_val),
    shuffle=True,
)

################################################################################
#! Import abnormal reads
################################################################################

df_ab = pd.read_csv(file_ab, header=None, sep="\t")
df_ab.columns = ["seqID", "seq", "barcodeID"]

if "MIDS" in file_ab:
    df_ab.seq = "MIDS=" + df_ab.seq
else:
    df_ab.seq = "ACGT=" + df_ab.seq

X_real = onehot_encode_seq(df_ab.seq)

################################################################################
#! Novelity (Anomaly) detection
################################################################################
# ===========================================================
# ? FC layer
# ===========================================================
model_ = Model(model.get_layer(index=0).input, model.get_layer(index=-2).output)
# model_.summary()
train_vector = model_.predict(X_train, verbose=0, batch_size=32)
predict_vector = model_.predict(X_real, verbose=0, batch_size=32)

#===========================================================
#? LocalOutlierFactor
#===========================================================

metrics =  ["cityblock", "cosine", "euclidean", "l1", "l2", "manhattan",
            "braycurtis", "canberra", "chebyshev", "correlation", "dice", "hamming", "jaccard", "kulsinski", "minkowski", "rogerstanimoto", "russellrao",
            "sokalmichener", "sokalsneath", "sqeuclidean", "yule"]

df_res = pd.DataFrame()
for metric in metrics:
    print("======================")
    print(metric)
    clf = LocalOutlierFactor(
    n_neighbors=20,
    metric=metric,
    leaf_size=30,
    novelty=True,
    n_jobs=threads,
    )
    clf.fit(train_vector)
    outliers = clf.predict(predict_vector)
    outliers = np.where(outliers == 1, "normal", "abnormal")
    df_ab["outliers"] = outliers
    print(df_ab.groupby("barcodeID").outliers.value_counts())
    df_test = df_ab.groupby("barcodeID").outliers.value_counts().to_frame()
    df_test["metric"] = metric
    df_res = pd.concat([df_res, df_test])

################################################################################
#! Evaluation
################################################################################

from sklearn.metrics import accuracy_score
# from sklearn.metrics import precision_score, recall_score, f1_score

df_ab["true"] = np.where(df_ab.barcodeID.str.contains("nega"), 1, 0)

df_report = pd.DataFrame(df_ab.groupby("barcodeID").apply(lambda x: accuracy_score(x.true, x.pred)))
df_report.columns = ["value"]
df_report["model"] = L2

if "MIDS" in file_ab:
    df_report["convertion"] = "MIDS"
else:
    df_report["convertion"] = "ACGT"

df_report.to_csv("accuracy_anomaly_detection.csv", mode='a', header=False)

