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
from tensorflow.keras import backend as K
from tensorflow.keras import regularizers, utils
from tensorflow.keras.regularizers import l2
from tensorflow.keras.layers import (Activation, BatchNormalization, Conv1D, Dense, Flatten,
                                    GlobalMaxPooling1D, MaxPooling1D,Dropout)
from tensorflow.keras.models import Model, Sequential, load_model
from tensorflow.keras.callbacks import EarlyStopping, ModelCheckpoint

os.environ['TF_KERAS'] = '1'
from keras_radam import RAdam

# ============================================================


###############################################
# Read real reads and one-hot encording
###############################################

args = sys.argv
file_name = args[1]

# file_name = '.DAJIN_temp/data/DAJIN_real.txt'

df_real = pd.read_csv(file_name, header=None, sep='\t')
df_real.columns = ["seqID", "seq", "barcodeID"]
df_real.seq = "MIDS=" + df_real.seq

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

print("One-hot encording to real reads...")
X_real = X_onehot(df_real.seq)

###############################################
# load model and X_test
###############################################

model = load_model(".DAJIN_temp/data/model_final.h5",
                   custom_objects={'RAdam': RAdam})

X_test = np.load(".DAJIN_temp/data/x_test.npz")
X_test = X_test["X_test"]

labels_index = pd.read_csv(".DAJIN_temp/data/labels_index.txt")
labels_index.columns = ["label"]
###############################################
# Read cos sim simulation file
###############################################

sim_score = pd.read_csv(".DAJIN_temp/data/cosine_sim.txt", header=None, names=["score"])
cossim_threshold = sim_score.score.quantile(0.0001)

###############################################
# Compute cosine similarity
###############################################

def get_score_cosine(model, reference, query):
    model_ = Model(model.get_layer(index=0).input,
                model.get_layer(index=-3).output)  # Delete FC layer
    # print(model_.summary())
    #print("Obtain L2-normalized vectors from the simulated reads...")
    normal_vector = model_.predict(reference, verbose=0, batch_size=64)
    #print("Obtain L2-normalized vectors of the real reads...")
    predict_vector = model_.predict(query, verbose=0, batch_size=64)
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

cos_all, normal_vector, predict_vector = get_score_cosine(
    model, X_test[:10000,], X_real)

# df_all = pd.concat([df_real.reset_index(
#    drop=True), pd.DataFrame(cos_all)], axis=1)
#df_all.columns = ["barcodeID", "seqID", "cos_similarity"]

df_anomaly = df_real[["barcodeID", "seqID"]]
df_anomaly["abnormal_prediction"] = pd.Series(cos_all).apply(
    lambda x: "normal" if x > cossim_threshold else "abnormal")


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

df_result.to_csv('.DAJIN_temp/data/DAJIN_prediction_result.txt',
                 sep='\t', index=False, header=False)
                 