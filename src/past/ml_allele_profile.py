import os
os.environ['TF_CPP_MIN_LOG_LEVEL']='3'
os.environ['TF_FORCE_GPU_ALLOW_GROWTH']='true'

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

import tensorflow as tf
from tensorflow.keras import regularizers, utils
from tensorflow.keras.layers import Activation, Conv1D, Dense, Flatten,MaxPooling1D
from tensorflow.keras.models import Model, Sequential, load_model
from tensorflow.keras.callbacks import EarlyStopping


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
# df_real = df_real[df_real.barcodeID=="barcode14"].reset_index(drop=True)

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
    # X = X[:, :, 1:]
    return(X)

print("One-hot encording simulated reads...")
X_sim = X_onehot(df_sim.seq)
print("One-hot encording real reads...")
X_real = X_onehot(df_real.seq)

# X_sim_reshape = X_sim.reshape([X_sim.shape[0], X_sim.shape[1], X_sim.shape[2], 1])
# X_real_reshape = X_real.reshape([X_real.shape[0], X_real.shape[1], X_real.shape[2], 1])
# X_real_reshape = X_real.reshape([X_real.shape[0], -1])
# ====================================
# ## Train test split
# ====================================
print("Model training...")

labels, labels_index = pd.factorize(df_sim.barcodeID)
labels_categorical = utils.to_categorical(labels)

X_train, X_test, Y_train, Y_test = train_test_split(
    X_sim, labels_categorical,
    test_size=0.2, shuffle=True)

# X_train, X_test, Y_train, Y_test = train_test_split(
#     X_sim_reshape.astype("float16"), labels_categorical,
#     test_size=0.2, shuffle=True)

# ====================================
# # Model construction
# ====================================

model = Sequential()

model.add(Conv1D(filters=32, kernel_size=32, activation="relu", padding = "same",
                input_shape=(X_train.shape[1], X_train.shape[2]), name="1st_Conv1D"))
model.add(MaxPooling1D(pool_size=4, name="1st_MaxPooling1D"))

model.add(Conv1D(filters=64, kernel_size=16, padding = "same",
                activation="relu", name="2nd_Conv1D"))
model.add(MaxPooling1D(pool_size=8, name="2nd_MaxPooling1D"))

model.add(Conv1D(filters=128, kernel_size=8, padding = "same",
                activation="relu", name="3rd_Conv1D"))
model.add(MaxPooling1D(pool_size=16, name="3rd_MaxPooling1D"))

model.add(Flatten(name="flatten"))

# model.add(Dense(32, activation='relu', name="1st_FC"))
# L2 <<<<<<<
alpha = 0.1  # hyperparameter
model.add(Dense(32, activation='linear',
                activity_regularizer=regularizers.l2(alpha), name="L2-normalization"))

model.add(Dense(len(labels_index), activation='softmax', name="final_layer"))
model.compile(optimizer='adam', loss='categorical_crossentropy',
            metrics=['accuracy'])
model.summary()
# tf.keras.utils.plot_model(model, show_shapes=True, show_layer_names=False, rankdir = "LR", to_file='model.png', dpi=350)
# -
early_stopping = EarlyStopping(monitor='val_loss', patience=10) 

stack = model.fit(X_train, Y_train, epochs=100, verbose=1, batch_size = 32,
                validation_split=0.2, shuffle=True,
                callbacks = [early_stopping])

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
    return score

def cosine_similarity(x1, x2):
    if x1.ndim == 1:
        x1 = x1[np.newaxis]
    if x2.ndim == 1:
        x2 = x2[np.newaxis]
    x1_norm = np.linalg.norm(x1, axis=1)
    x2_norm = np.linalg.norm(x2, axis=1)
    cosine_sim = np.dot(x1, x2.T)/(x1_norm*x2_norm+1e-10)
    return cosine_sim

print("Abnormal allele detection...")
cos_all = get_score_cosine(model, X_train[0:500], np.concatenate([X_sim, X_real]))
# del X_sim

df_anomaly = pd.concat([df_sim[["barcodeID", "seqID"]].reset_index(drop=True),
                    df_real[["barcodeID", "seqID"]].reset_index(drop=True)])

df_anomaly["cos_sim"] = cos_all

df_anomaly["label"] = df_anomaly.barcodeID.apply(
    lambda x: 1 if "simulated" in x else -1)

optimal_threshold = df_anomaly[df_anomaly.label == 1].cos_sim.quantile(0.001)

df_anomaly["abnormal_prediction"] = (
    df_anomaly[df_anomaly.label == -1].reset_index().
    cos_sim.apply(lambda x: "normal" if x > optimal_threshold else "abnormal")
)
df_anomaly.drop(["cos_sim", "label"], axis=1, inplace=True)
df_anomaly = df_anomaly[~df_anomaly.barcodeID.str.contains("simulated")].reset_index(drop=True)
# df_anomaly.groupby("barcodeID").abnormal_prediction.value_counts()
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

# del X_real

df_predict = pd.Series(predict, dtype="str") + "_"

for i, j in enumerate(labels_index.str.replace("_simulated.*$", "")):
    df_predict = df_predict.str.replace(str(i)+"_", j)

df_result = pd.DataFrame({"barcodeID": df_anomaly.iloc[:, 0],
                          "seqID": df_anomaly.iloc[:, 1],
                          "predict": df_predict,
                          "anomaly": df_anomaly.iloc[:, 2]})

df_result.predict = df_result.predict.mask(
    df_result.anomaly.str.contains("abnormal"), df_result.anomaly)

df_result = df_result.drop('anomaly', axis=1)
# df_result.predict.value_counts()

# ====================================
# Output the results
# ====================================

df_result.to_csv('.DAJIN_temp/data/DAJIN_MIDS_prediction_result.txt',
                 sep='\t', index=False, header=False)

