import os

os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"
os.environ["TF_FORCE_GPU_ALLOW_GROWTH"] = "true"

import sys
import random as rn
import numpy as np
import pandas as pd

import pickle

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder, LabelBinarizer
from sklearn.neighbors import LocalOutlierFactor

import tensorflow as tf
from tensorflow.keras.layers import Input, Conv1D, Dense, Flatten, MaxPooling1D
from tensorflow.keras.models import Model

# set random seed
os.environ['PYTHONHASHSEED'] = '0'
np.random.seed(111)
rn.seed(12345)
tf.random.set_seed(1234)

################################################################################
#! I/O naming
################################################################################

#===========================================================
#? Auguments
#===========================================================

args = sys.argv
file_name = args[1]
threads = int(args[2])

if threads == "":
    import multiprocessing
    threads = multiprocessing.cpu_count() // 2

#===========================================================
#? Input
#===========================================================

df_sim = pd.read_csv(file_name, header=None, sep="\t")
df_sim.columns = ["seqID", "seq", "barcodeID"]
df_sim.seq = "MIDS=" + df_sim.seq

################################################################################
#! Training model
################################################################################

#===========================================================
#? Define Function
#===========================================================


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


#===========================================================
#? Train test split
#===========================================================

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

#===========================================================
#? Training
#===========================================================

callback = tf.keras.callbacks.EarlyStopping(monitor='val_loss', patience=20)

model.fit(
    X_train,
    Y_train,
    epochs=200,
    verbose=0,
    batch_size=32,
    validation_data=(X_val, Y_val),
    callbacks=[callback],
    shuffle=True,
)

model_eval = model.evaluate(X_val, Y_val, verbose=0)

print(f"model loss: {model_eval[0]}", file=sys.stderr)
print(f"model accuracy: {model_eval[1]}", file=sys.stderr)

################################################################################
#! Novelty (Anomaly) detection
################################################################################
#===========================================================
#? Extract layer
#===========================================================

model_ = Model(model.get_layer(index=0).input, model.get_layer(index=-2).output)
train_vector = model_.predict(X_train, verbose=0, batch_size=32)

#===========================================================
#? LocalOutlierFactor
#===========================================================

clf = LocalOutlierFactor(
    n_neighbors=100,
    metric="euclidean",
    contamination="auto",
    leaf_size=30,
    novelty=True,
    n_jobs=threads,
)

clf.fit(train_vector)

################################################################################
#! Save models
################################################################################


print(f"loss: {model_eval[0]}", file=open('.DAJIN_temp/data/model_evaluate.csv', 'a'))
print(f"accuracy: {model_eval[1]}", file=open('.DAJIN_temp/data/model_evaluate.csv', 'a'))

np.save(".DAJIN_temp/data/labels_index.npy", labels_index)

model.save(".DAJIN_temp/data/model.h5")

pickle.dump(clf, open(".DAJIN_temp/data/model_lof.sav", "wb"))

