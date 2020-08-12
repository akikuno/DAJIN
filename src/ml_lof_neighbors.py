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

# file_train = ".DAJIN_temp/data/DAJIN_MIDS_train.txt"
# file_test = ".DAJIN_temp/data/DAJIN_MIDS_test.txt"
# threads = 12

#==========================================================
#? Auguments
#==========================================================

args = sys.argv
file_train = args[1]
file_test = args[2]
if args[3] == "":
    import multiprocessing
    threads = multiprocessing.cpu_count() // 2
else:
    threads = int(args[3])

iter_num = int(args[4])

print(f"{iter_num} execution...")

#==========================================================
#? Input
#==========================================================

df_sim = pd.read_csv(file_train, header=None, sep="\t")
df_sim.columns = ["seqID", "seq", "barcodeID"]

if "MIDS" in file_train:
    df_sim.seq = "MIDS=" + df_sim.seq
else:
    df_sim.seq = "ACGT=" + df_sim.seq

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

df_test = pd.read_csv(file_test, header=None, sep="\t")
df_test.columns = ["seqID", "seq", "barcodeID"]

if "MIDS" in file_test:
    df_test.seq = "MIDS=" + df_test.seq
else:
    df_test.seq = "ACGT=" + df_test.seq

X_real = onehot_encode_seq(df_test.seq)

################################################################################
#! Novelity (Anomaly) detection
################################################################################

# ===========================================================
# ? FC layer
# ===========================================================

model_ = Model(model.get_layer(index=0).input, model.get_layer(index=-2).output)
train_vector = model_.predict(X_train, verbose=0, batch_size=32)
predict_vector = model_.predict(X_real, verbose=0, batch_size=32)

#===========================================================
#? LocalOutlierFactor
#===========================================================

df_tmp = df_test

metrics =  ["cityblock", "cosine", "euclidean", "l1", "l2", "manhattan",
            "braycurtis", "canberra", "chebyshev", "correlation", "dice", "hamming", "jaccard", "kulsinski", "minkowski", "rogerstanimoto", "russellrao",
            "sokalmichener", "sokalsneath", "sqeuclidean"]

# metrics = ["cosine", "jaccard"]
df_res = pd.DataFrame()
for n_neighbor in [5, 10, 20, 50, 100, 200, 500, 1000]:
    clf = LocalOutlierFactor(
        n_neighbors=n_neighbor,
        metric="euclidean",
        novelty=True,
        n_jobs=threads,
    )
    clf.fit(train_vector)
    outliers = clf.predict(predict_vector)
    outliers = np.where(outliers == 1, "normal", "abnormal")
    df_test["outliers"] = outliers
    df_tmp = df_test.groupby("barcodeID").outliers.value_counts().to_frame()
    df_tmp["n_neighbor"] = n_neighbor
    df_tmp["iter"] = iter_num
    df_res = pd.concat([df_res, df_tmp])

df_res.to_csv("results_lof_neighbor.csv", mode = "a", header=False)