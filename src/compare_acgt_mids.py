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

import umap

################################################################################
#! I/O naming
################################################################################

#==========================================================
#? Arguments
#==========================================================
key = "MIDS"

file_train = f".DAJIN_temp/data/DAJIN_{key}_train.txt"
file_test = f".DAJIN_temp/data/DAJIN_{key}_test.txt"

df_train = pd.read_csv(file_train, header=None, sep="\t")
df_test = pd.read_csv(file_test, header=None, sep="\t")

for i in [df_train, df_test]:
    i.columns = ["seqID", "seq", "barcodeID"]

df_train.seq += f"{key}="
df_test.seq += f"{key}="


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

labels, labels_index = pd.factorize(df_train.barcodeID)
labels = tf.keras.utils.to_categorical(labels)

X_train, X_val, Y_train, Y_val = train_test_split(
    onehot_encode_seq(df_train.seq), labels, test_size=0.2, shuffle=True
)


#==========================================================
#? CNN
#==========================================================

inputs = Input(shape = (X_train.shape[1], X_train.shape[2]))
init_kernel_size = int(128)

x = Conv1D(
        filters=32,
        kernel_size=init_kernel_size,
        activation="relu",
        name="1st_Conv1D",
    )(inputs)
x = MaxPooling1D(pool_size=4, name="1st_MaxPooling1D")(x)

x = Conv1D(
        filters=32,
        kernel_size=int(init_kernel_size / 2),
        activation="relu",
        name="2nd_Conv1D",
    )(x)
x = MaxPooling1D(pool_size=4, name="2nd_MaxPooling1D")(x)

x = Conv1D(
        filters=32,
        kernel_size=int(init_kernel_size / 4),
        activation="relu",
        name="3rd_Conv1D",
    )(x)
x = MaxPooling1D(pool_size=4, name="3rd_MaxPooling1D")(x)

x = Conv1D(
        filters=32,
        kernel_size=int(init_kernel_size / 16),
        activation="relu",
        name="4th_Conv1D",
    )(x)
x = MaxPooling1D(pool_size=4, name="4th_MaxPooling1D")(x)

x = Flatten(name="flatten")(x)

x = Dense(64, activation="relu", name="1st_FC")(x)

# if L2 == "on":
#     alpha = 0.1
#     x = alpha * tf.divide(x, tf.norm(x, ord='euclidean'))
# else:
#     pass

predictions = Dense(len(labels_index), activation="softmax", name="softmax")(x)

model = Model(inputs = inputs, outputs = predictions)

model.compile(optimizer="adam", loss="categorical_crossentropy", metrics=["accuracy"])
model.summary()

#==========================================================
#? Training
#==========================================================

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
#! Novelity (Anomaly) detection
################################################################################

#==========================================================
#? last layer
#==========================================================

model_ = Model(model.get_layer(index=0).input, model.get_layer(index=-2).output)
train_vector = model_.predict(X_train, verbose=0, batch_size=32)

# del X_train  # <<<

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
#! Novelity (Anomaly) detection
################################################################################

# ===========================================================
# ? L2 layer
# ===========================================================
X_test = onehot_encode_seq(df_test.seq)

predict_vector = model_.predict(X_test, verbose=0, batch_size=32)

# ===========================================================
# ? LocalOutlierFactor
# ===========================================================

outliers = clf.predict(predict_vector)
outliers = np.where(outliers == 1, 1, 0)

df_test["pred"] = outliers

df_test.groupby("barcodeID").pred.value_counts()

################################################################################
#! UMAP
################################################################################

#===========================================================
#? embedding
#===========================================================

umap_before = umap.UMAP().fit_transform(X_test.reshape(X_test.shape[0],-1))
umap_after = umap.UMAP().fit_transform(predict_vector)

df_umap_before = pd.DataFrame(umap_before)
df_umap_after = pd.DataFrame(umap_after)

df_umap_before["id"] = "before"
df_umap_after["id"] = "after"
df_umap_before["label"] = df_test.barcodeID
df_umap_after["label"] = df_test.barcodeID
df_umap_before["convertion"] = key
df_umap_after["convertion"] = key

df_umap = pd.concat([df_umap_test_before, df_umap_test_after], axis=0)


df_umap.to_csv(f"{key}_umap.csv")
