import os

os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"
os.environ["TF_FORCE_GPU_ALLOW_GROWTH"] = "true"

import sys
import numpy as np
import pandas as pd

import pickle

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder, LabelBinarizer
from sklearn.neighbors import LocalOutlierFactor

import tensorflow as tf
from tensorflow.keras import regularizers
from tensorflow.keras.layers import Conv1D, Dense, Flatten, MaxPooling1D
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

L2 = args[4]

print(args)

if L2 == "on" or L2 == "off":
    pass
else:
    raise ValueError("Invalid parameter in L2")

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

X_train, X_test, Y_train, Y_test = train_test_split(
    onehot_encode_seq(df_sim.seq), labels, test_size=0.2, shuffle=True
)

del df_sim
del X_test
del Y_test

#==========================================================
#? L2-constrained Softmax Loss
#==========================================================

init_kernel_size = int(128)

model = tf.keras.Sequential()
model.add(
    Conv1D(
        filters=32,
        kernel_size=init_kernel_size,
        activation="relu",
        input_shape=(X_train.shape[1], X_train.shape[2]),
        name="1st_Conv1D",
    )
)
model.add(MaxPooling1D(pool_size=4, name="1st_MaxPooling1D"))

model.add(
    Conv1D(
        filters=32,
        kernel_size=int(init_kernel_size / 2),
        activation="relu",
        name="2nd_Conv1D",
    )
)
model.add(MaxPooling1D(pool_size=4, name="2nd_MaxPooling1D"))

model.add(
    Conv1D(
        filters=32,
        kernel_size=int(init_kernel_size / 4),
        activation="relu",
        name="3rd_Conv1D",
    )
)
model.add(MaxPooling1D(pool_size=4, name="3rd_MaxPooling1D"))

model.add(
    Conv1D(
        filters=32,
        kernel_size=int(init_kernel_size / 16),
        activation="relu",
        name="4th_Conv1D",
    )
)
model.add(MaxPooling1D(pool_size=4, name="4th_MaxPooling1D"))

model.add(Flatten(name="flatten"))

# model.add(Dense(64, activation="relu", name="1st_FC"))

if L2 == "on":
    alpha = 0.1
    model.add(
        Dense(
            32,
            activation="linear",
            activity_regularizer=regularizers.l2(alpha),
            name="L2",
        )
    )
else:
    model.add(
        Dense(
            32,
            activation="linear",
            name="linear",
        )
    )

model.add(Dense(len(labels_index), activation="softmax", name="softmax"))
model.compile(optimizer="adam", loss="categorical_crossentropy", metrics=["accuracy"])
model.summary()

#==========================================================
#? Training
#==========================================================

model.fit(
    X_train,
    Y_train,
    epochs=20,
    verbose=0,
    batch_size=32,
    validation_split=0.2,
    shuffle=True,
)

################################################################################
#! Novelity (Anomaly) detection
################################################################################

#==========================================================
#? L2 layer
#==========================================================

model_ = Model(model.get_layer(index=0).input, model.get_layer(index=-2).output)
train_vector = model_.predict(X_train, verbose=0, batch_size=32)

del X_train  # <<<

#==========================================================
#? LocalOutlierFactor
#==========================================================

clf = LocalOutlierFactor(
    n_neighbors=20,
    metric="euclidean",
    contamination="auto",
    leaf_size=400,
    novelty=True,
    n_jobs=threads,
)

clf.fit(train_vector)


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
# ? L2 layer
# ===========================================================

predict_vector = model_.predict(X_real, verbose=0, batch_size=32)

# ===========================================================
# ? LocalOutlierFactor
# ===========================================================

outliers = clf.predict(predict_vector)
outliers = np.where(outliers == 1, 1, 0)

df_ab["pred"] = outliers

df_ab.groupby("barcodeID").pred.value_counts()

df_ab["true"] = np.where(df_ab.barcodeID.str.contains("nega"), 1, 0)

################################################################################
#! Evaluation
################################################################################

from sklearn.metrics import accuracy_score
# from sklearn.metrics import precision_score, recall_score, f1_score

df_report = pd.DataFrame(df_ab.groupby("barcodeID").apply(lambda x: accuracy_score(x.true, x.pred)))
df_report.columns = ["value"]
df_report["model"] = L2

if "MIDS" in file_ab:
    df_report["convertion"] = "MIDS"
else:
    df_report["convertion"] = "ACGT"

df_report.to_csv("accuracy_anomaly_detection.csv", mode='a', header=False)
# print(df_ab.groupby("barcodeID").apply(lambda x: precision_score(x.true, x.pred)))
# print(df_ab.groupby("barcodeID").apply(lambda x: recall_score(x.true, x.pred)))
# print(df_ab.groupby("barcodeID").apply(lambda x: f1_score(x.true, x.pred)))

################################################################################
#! Save models
################################################################################
# np.save(".DAJIN_temp/data/labels_index.npy", labels_index)

# model.save(".DAJIN_temp/data/model_l2.h5")

# pickle.dump(clf, open(".DAJIN_temp/data/model_lof.sav", "wb"))

