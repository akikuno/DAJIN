import os
os.environ['TF_CPP_MIN_LOG_LEVEL']='3'
os.environ['TF_KERAS'] = '1'
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

import scipy.sparse

import tensorflow as tf
from tensorflow.keras.optimizers import Adam, SGD
from tensorflow.keras import backend as K
from tensorflow.keras import regularizers, utils
from tensorflow.keras.regularizers import l2
from tensorflow.keras.layers import (Activation, BatchNormalization, Conv1D, Dense, Flatten, LSTM, 
                                    GlobalMaxPooling1D, MaxPooling1D,Dropout, Embedding)
from tensorflow.keras.models import Model, Sequential, load_model

from keras_radam import RAdam
from tensorflow.keras.callbacks import EarlyStopping 


###############################################
# ハイバーパラメータ
###############################################
batch_size = 64

###############################################
# モデルの構築
###############################################
model = Sequential()
model.add(Embedding(20000, 128))
model.add(LSTM(128, dropout=0.2, recurrent_dropout=0.2))
model.add(Dense(6, activation='sigmoid'))

model.compile(loss='categorical_crossentropy',
              optimizer='adam',
              metrics=['accuracy'])
###############################################
# ジェネレータとステップ数
###############################################
from typing import List

def process_input(input_data):
    # data増強をこちらに
    return input_data

def X_onehot(X_data):
    sequence = X_data[i]
    seq_array = np.array(list(sequence))
    label_encoder = LabelEncoder()
    integer_encoded_seq = label_encoder.fit_transform(seq_array)
    # one hot the sequence
    onehot_encoder = OneHotEncoder(sparse=False, categories='auto', dtype='uint8')
    integer_encoded_seq = integer_encoded_seq.reshape(
        len(integer_encoded_seq), 1)
    onehot_encoded_seq = onehot_encoder.fit_transform(integer_encoded_seq)
    X = onehot_encoded_seq[:, :, 1:]
    return(X)

def count_lines(file_path:str):
    count=0
    with open(file_path, 'r') as f:
        for line in f:
            count += 1
    return count

# ここでは、ラベル作成のために必要であろうと思ってclass_mapも追加しておきました、ラベルがすでに数字がされていましたら要りませんです！
def xtrain_generator(path: str, batch_size=1000):
        i = 0
        max_i = count_lines(path)
        while True:
            # リセットします
            X_sim = np.empty((0, 2729, 5),
                        dtype="uint8")
            labels_df = pd.DataFrame()
            add_n = batch_size
            # たまに、うまく行かないのもあるので、それを外て埋め直したいために、whileに入れる
            while len(X_sim) < batch_size:
                # batch_sizeに合わせた行数を読みます。skiprows=range(1, i)は項目名(最初のrow)を読むためです.
                df_sim = pd.read_csv(path,
                                        dtype=str,
                                        skiprows=range(1, i),
                                        nrows=add_n)
                df_sim.columns = ["seqID", "seq", "barcodeID"]
                df_sim.seq = "MIDS=" + df_sim.seq
                # -------------------------------------------------
                # # input: One-hot encording
                # -------------------------------------------------
                X_sim = np.append(X_sim, X_onehot(df_sim.seq), axis=0)
                # -------------------------------------------------
                # # label
                # -------------------------------------------------
                labels, labels_index = pd.factorize(df_sim.barcodeID)
                labels_categorical = utils.to_categorical(labels)
                labels_df = pd.concat([labels_df, labels_categorical])
                i += add_n
                if i > max_i:  # ファイルが終わったら最初からやり直す
                    i = 0
                add_n = batch_size - len(inputs)
            # -------------------------------------------------
            # # Train test split
            # -------------------------------------------------
            X_train, X_test, Y_train, Y_test = train_test_split(
                X_sim, labels_df,
                test_size=0.2, shuffle=True)
            yield X_train


def ytrain_generator(path: str, batch_size=1000):
        i = 0
        max_i = count_lines(path)
        while True:
            # リセットします
            X_sim = np.empty((0, 2729, 5),
                        dtype="uint8")
            labels_df = pd.DataFrame()
            add_n = batch_size
            # たまに、うまく行かないのもあるので、それを外て埋め直したいために、whileに入れる
            while len(X_sim) < batch_size:
                # batch_sizeに合わせた行数を読みます。skiprows=range(1, i)は項目名(最初のrow)を読むためです.
                df_sim = pd.read_csv(path,
                                        dtype=str,
                                        skiprows=range(1, i),
                                        nrows=add_n)
                df_sim.columns = ["seqID", "seq", "barcodeID"]
                df_sim.seq = "MIDS=" + df_sim.seq
                # -------------------------------------------------
                # # input: One-hot encording
                # -------------------------------------------------
                X_sim = np.append(X_sim, X_onehot(df_sim.seq), axis=0)
                # -------------------------------------------------
                # # label
                # -------------------------------------------------
                labels, labels_index = pd.factorize(df_sim.barcodeID)
                labels_categorical = utils.to_categorical(labels)
                labels_df = pd.concat([labels_df, labels_categorical])
                i += add_n
                if i > max_i:  # ファイルが終わったら最初からやり直す
                    i = 0
                add_n = batch_size - len(inputs)
            # -------------------------------------------------
            # # Train test split
            # -------------------------------------------------
            X_train, X_test, Y_train, Y_test = train_test_split(
                X_sim, labels_df,
                test_size=0.2, shuffle=True)
            yield Y_train

# 大量のデータで学習しているからこそチェックポイントなどをつけましょう。
ckpt = tf.keras.callbacks.ModelCheckpoint('.DAJIN_temp/data/ckpt.hdf5', monitor='val_loss',
                                          verbose=0, save_best_only=True,
                                          save_weights_only=False, mode='auto', period=3)
lr = tf.keras.callbacks.ReduceLROnPlateau(monitor='val_loss', factor=0.1, patience=5,
                                          verbose=1, mode='auto', min_delta=0.0001,
                                          cooldown=3, min_lr=0)

early_stopping = EarlyStopping(monitor='val_loss', patience=10)

model.fit_generator(
    xtrain_generator("DAJIN_test_sim.txt"),
    steps_per_epoch=800 // batch_size,
    validation_data=ytrain_generator(path="DAJIN_test_sim.txt"),
    validation_steps= 200 // batch_size,
    epochs=30,           
    callbacks=[ckpt, lr, early_stopping]
)