import os
os.environ['TF_CPP_MIN_LOG_LEVEL']='3'
import warnings
warnings.filterwarnings("ignore")
warnings.simplefilter(action='ignore', category=FutureWarning)

import sys
import re
from functools import partial
from tqdm import tqdm

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import OneHotEncoder

from tensorflow.keras import backend as K
from tensorflow.keras import regularizers, utils
from tensorflow.keras.layers import (Activation, Conv1D, Dense, Flatten,
                                    MaxPooling1D)
from tensorflow.keras.models import Model, Sequential, load_model


# ====================================
# Input and format data
# ====================================

args = sys.argv
file_name = args[1]
# file_name = ".DAJIN_temp/data/DAJIN_MIDS.txt"

df = pd.read_csv(file_name, header=None, sep='\t')
df.columns = ["seqID", "seq", "barcodeID"]
df.seq = "MIDS=" + df.seq

df_sim = df[df.barcodeID.str.contains("simulated")].reset_index(drop=True)
df_real = df[~df.barcodeID.str.contains("simulated")].reset_index(drop=True)

# Output names
fig_dirs = ["results/figures/png", "results/figures/svg"]
# output_npz = file_name.replace(".txt", ".npz")
# output_figure = file_name.replace(".txt", "").replace("data_for_ml/", "")
# output_model = file_name.replace(".txt", ".h5")

# # ====================================
# # One-hot encording
# # ====================================
def X_onehot(X_data):
    X = np.empty((len(X_data), len(X_data[0]), 5),
        dtype="uint8")
    for i in tqdm(range(0, len(X_data))):
        sequence = X_data[i]
        seq_array = np.array(list(sequence))
        label_encoder = LabelEncoder()
        integer_encoded_seq = label_encoder.fit_transform(seq_array)
        # one hot the sequence
        onehot_encoder = OneHotEncoder(sparse=False, categories='auto', dtype='uint8')
        integer_encoded_seq = integer_encoded_seq.reshape(
            len(integer_encoded_seq), 1)
        onehot_encoded_seq = onehot_encoder.fit_transform(integer_encoded_seq)
        X[i] = onehot_encoded_seq
    X = X[:, :, 1:]
    return(X)

print("One-hot encording simulated reads...")
X_sim = X_onehot(df_sim.seq)
print("One-hot encording real reads...")
X_real = X_onehot(df_real.seq)

# ====================================
# ## Train test split
# ====================================
print("Model training...")

labels, labels_index = pd.factorize(df_sim.barcodeID)
labels_categorical = utils.to_categorical(labels)

X_train, X_test, Y_train, Y_test = train_test_split(
    X_sim, labels_categorical,
    test_size=0.2, shuffle=True)

# ====================================
# # Model construction
# ====================================

model = Sequential()

model.add(Conv1D(filters=32, kernel_size=32, activation="relu",
                input_shape=(X_sim.shape[1], X_sim.shape[2]), name="1st_Conv1D"))
model.add(MaxPooling1D(pool_size=4, name="1st_MaxPooling1D"))

model.add(Conv1D(filters=32, kernel_size=16,
                activation="relu", name="2nd_Conv1D"))
model.add(MaxPooling1D(pool_size=4, name="2nd_MaxPooling1D"))

model.add(Conv1D(filters=32, kernel_size=8,
                activation="relu", name="3rd_Conv1D"))
model.add(MaxPooling1D(pool_size=4, name="3rd_MaxPooling1D"))

model.add(Conv1D(filters=32, kernel_size=4,
                activation="relu", name="4th_Conv1D"))
model.add(MaxPooling1D(pool_size=4, name="4th_MaxPooling1D"))

model.add(Flatten(name="flatten"))

model.add(Dense(64, activation='relu', name="1st_FC"))

# L2 <<<<<<<
alpha = 0.1  # hyperparameter
model.add(Dense(64, activation='linear',
                activity_regularizer=regularizers.l2(alpha), name="L2-softmax"))

model.add(Dense(len(labels_index), activation='softmax', name="final_layer"))
model.compile(optimizer='adam', loss='categorical_crossentropy',
            metrics=['accuracy'])
# model.summary()
# -
stack = model.fit(X_train, Y_train, epochs=20, verbose=0,
                validation_split=0.2, shuffle=True)


# evaluate = model.evaluate(X_test, Y_test, verbose=0)


# ====================================
# ## Compute cosine similarity
# ====================================

def get_score_cosine(model, train, test):
    model_ = Model(model.get_layer(index=0).input,
                model.get_layer(index=-2).output)  # Delete FC layer
    # print(model_.summary())
    print("Obtain L2-normalized vectors from the simulated reads...")
    normal_vector = model_.predict(train, verbose=1, batch_size=32)
    print("Obtain L2-normalized vectors of the real reads...")
    predict_vector = model_.predict(test, verbose=1, batch_size=32)
    score_len = len(predict_vector)
    score = np.zeros(score_len)
    print("Calculate cosine similarity to detect abnormal allele ...")
    for i in tqdm(range(score_len)):
        score[i] = cosine_similarity(predict_vector[i], normal_vector).max()
    return score, normal_vector, predict_vector

def cosine_similarity(x1, x2):
    if x1.ndim == 1:
        x1 = x1[np.newaxis]
    if x2.ndim == 1:
        x2 = x2[np.newaxis]
    x1_norm = np.linalg.norm(x1, axis=1)
    x2_norm = np.linalg.norm(x2, axis=1)
    cosine_sim = np.dot(x1, x2.T)/(x1_norm*x2_norm+1e-10)
    return cosine_sim


X_all = np.concatenate([X_sim, X_real])
print("Abnormal allele detection...")
cos_all, normal_vector, predict_vector = get_score_cosine(
    model, X_train[0:1000], X_all)

df_name = pd.concat([df_sim[["barcodeID", "seqID"]].reset_index(drop=True),
                    df_real[["barcodeID", "seqID"]].reset_index(drop=True)])

df_all = pd.concat([df_name.reset_index(
    drop=True), pd.DataFrame(cos_all)], axis=1)
df_all.columns = ["barcodeID", "seqID", "cos_similarity"]

df_all["label"] = df_all.barcodeID.apply(
    lambda x: 1 if "simulated" in x else -1)

optimal_threshold = df_all[df_all.label == 1].cos_similarity.quantile(0.001)

df_all["abnormal_prediction"] = df_all.cos_similarity.apply(
    lambda x: 1 if x > optimal_threshold else - 1)
    
df_anomaly = df_real[["barcodeID", "seqID"]]
df_anomaly["abnormal_prediction"] = df_all[df_all.label == -1].reset_index().cos_similarity.apply(
    lambda x: "normal" if x > optimal_threshold else "abnormal")

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
# Output the results
# ====================================

df_result.to_csv('.DAJIN_temp/data/DAJIN_prediction_result.txt',
                 sep='\t', index=False, header=False)