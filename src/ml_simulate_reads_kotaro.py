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
from tensorflow.keras.optimizers import Adam
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

args = sys.argv
file_name = args[1]
# file_name = 'drive/My Drive/DAJIN_materials/pointmutation/test.txt.gz'
# file_name = 'DAJIN_sim_test.txt'
# file_name = 'DAJIN_sim.txt'

df_sim = pd.read_csv(file_name, header=None, sep='\t')
df_sim.columns = ["seqID", "seq", "barcodeID"]
df_sim.seq = "MIDS=" + df_sim.seq

# Output names
# output_name = file_name.replace(".txt", "")

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

print("One-hot encording to simulated reads...")
X_sim = X_onehot(df_sim.seq)

# np.savez_compressed('.DAJIN_temp/data/x_sim.npz', X_sim=X_sim)
# X_sim = np.load('.DAJIN_temp/data/x_sim.npz')
# X_sim = X_sim["X_sim"]

###############################################
# Train-test split
###############################################

labels, labels_index = pd.factorize(df_sim.barcodeID)
pd.Series(labels_index.values).to_csv(".DAJIN_temp/data/labels_index.txt", index=False)

labels_categorical = utils.to_categorical(labels)

X_train, X_test, Y_train, Y_test = train_test_split(
    X_sim, labels_categorical,
    test_size=0.2, shuffle=True)
del X_sim

###############################################
# Model constraction
###############################################
model = Sequential()

model.add(Conv1D(filters=32, kernel_size=4, padding="same",
                input_shape=(X_train.shape[1], X_train.shape[2]),
                 name="1st_Conv1D",
                kernel_regularizer=l2(0.001)))
model.add(BatchNormalization())
model.add(Activation('relu'))
model.add(MaxPooling1D(pool_size=4, name="1st_MaxPooling1D"))

model.add(Conv1D(filters=64, kernel_size=4, padding = "same",
                name="2nd_Conv1D",kernel_regularizer=l2(0.01)))
model.add(BatchNormalization())
model.add(Activation('relu'))
model.add(MaxPooling1D(pool_size=2, name="2nd_MaxPooling1D"))

model.add(Flatten(name="flatten"))

model.add(Dense(128, activation='relu'))

alpha = 0.0000001  # hyperparameter
alpha = 0.0001  # hyperparameter
model.add(Dense(128, activation='linear',
                activity_regularizer=regularizers.l2(alpha), name="L2-normalization"))
model.add(Dense(len(labels_index), activation='softmax', name="final_layer",
               kernel_regularizer=l2(0.001),
                 bias_regularizer=l2(0.001)))

# optimizer = Adam()
optimizer = RAdam()
model.compile(loss='categorical_crossentropy',
              optimizer=optimizer,
              metrics=['accuracy'])

# model.summary()
# -
# if(os.path.exists('.DAJIN_temp/data/model_callback.h5')):
#     model.load_weights(".DAJIN_temp/data/model_callback.h5")

early_stopping = EarlyStopping(monitor='val_loss', patience=10) 
# checkpoint = ModelCheckpoint(
#     filepath=".DAJIN_temp/data/model_callback.h5",
#     save_best_only=True)

###############################################
# Training
###############################################
history = model.fit(X_train, Y_train, epochs=100, verbose=1, batch_size = 64,
                validation_split=0.2, shuffle=True,
                callbacks = [early_stopping])

model.save('.DAJIN_temp/data/model_final.h5')

test_loss, test_acc = model.evaluate(X_test, Y_test, verbose=0)
print(f'test loss: {test_loss}, test acc: {test_acc}')

###############################################
# Cosine similarity
###############################################

def get_score_cosine(model, train, test):
    model_ = Model(model.get_layer(index=0).input,
                model.get_layer(index=-2).output)  # Delete FC layer
    # print(model_.summary())
    #print("Obtain L2-normalized vectors from the simulated reads...")
    normal_vector = model_.predict(train, verbose=0, batch_size=64)
    #print("Obtain L2-normalized vectors of the real reads...")
    predict_vector = model_.predict(test, verbose=0, batch_size=64)
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

# print("Abnormal allele detection...")
cos_all, normal_vector, predict_vector = get_score_cosine(
    model, X_test[0: 10000, ], X_train[0: 10000, ])
# pd.Series(cos_all).describe()
np.savetxt('.DAJIN_temp/data/cosine_sim.txt', cos_all)

# with open('.DAJIN_temp/data/cosine_sim.txt', 'a') as file:
#     np.savetxt(file, cos_all)

np.savez_compressed('.DAJIN_temp/data/x_test.npz', X_test=X_test[0:10000,])

# # Confusion matrix
# from sklearn.metrics import confusion_matrix, accuracy_score
# Y_pred = model.predict(X_test)
# Y_pred = np.argmax(Y_pred, axis=1)
# cm = confusion_matrix(Y_test.argmax(1), Y_pred, normalize = "true")

# # Transform to df for easier plotting
# cm_df = pd.DataFrame(cm,
#                      index = labels_index, 
#                      columns = labels_index)

# plt.style.use('ggplot')
# fig, ax = plt.subplots(3,1,  figsize=(5, 15), facecolor='w')
# plt.tight_layout()

# fig.suptitle(f'RAdam', fontsize=20)
# fig.subplots_adjust(top=0.9)

# # summarize history for accuracy
# ax[0].plot(history.history['accuracy'])
# ax[0].plot(history.history['val_accuracy'])
# ax[0].set_title('Accuracy')
# ax[0].set_ylabel('accuracy')
# ax[0].set_xlabel('epoch')
# ax[0].legend(['train', 'validation'], loc='lower right')
# # summarize history for loss
# ax[1].plot(history.history['loss'])
# ax[1].plot(history.history['val_loss'])
# ax[1].set_title('Loss')
# ax[1].set_ylabel('loss')
# ax[1].set_xlabel('epoch')
# ax[1].legend(['train', 'validation'], loc='upper right')

# sns.heatmap(cm_df, annot=True, ax=ax[2])
# ax[2].set_ylabel('True label')
# ax[2].set_xlabel('Predicted label')
# plt.show()

# plt.savefig(f'test.jpg')


