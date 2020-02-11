from keras.models import Model
from keras import Model
from keras.models import load_model
import keras
from matplotlib.axes._axes import _log as matplotlib_axes_logger
import sys
import itertools
from sklearn.metrics import confusion_matrix
from keras.models import Sequential
from keras import regularizers
from keras.layers import Conv1D, Dense, MaxPooling1D, Flatten, Dropout, Activation
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from keras.utils import np_utils
from sklearn.model_selection import train_test_split
from functools import partial
import os
import numpy as np
import pandas as pd
import re
from sklearn.preprocessing import OneHotEncoder
from sklearn.preprocessing import LabelEncoder
from numpy import argmax
from numpy import array
from keras.backend import tensorflow_backend
from tensorflow.python.client import device_lib
import tensorflow as tf
import warnings
from tqdm import tqdm
warnings.filterwarnings('ignore')

# ====================================
# Input and format data
# ====================================

args = sys.argv
file_name = args[1]
# file_name = "data_for_ml/sequence_MIDS.txt.gz"
fig_dirs = ["results/figures/png", "results/figures/svg"]

for fig_dir in fig_dirs:
    os.makedirs(fig_dir, exist_ok=True)

os.makedirs("data_for_ml/model", exist_ok=True)
os.makedirs("results", exist_ok=True)

output_npz = file_name.replace(".txt.gz", ".npz").replace(
    "data_for_ml/", "data_for_ml/model/")
output_figure = file_name.replace(".txt.gz", "").replace("data_for_ml/", "")
output_model = file_name.replace(".txt.gz", "").replace(
    "data_for_ml", "data_for_ml/model")

df = pd.read_csv(file_name, header=None, sep='\t')
df.columns = ["barcodeID", "sequence", "seqID"]

if "ACGT" in file_name:
    df.sequence = "ACGT=" + df.sequence
    # print("ACGT")
else:
    df.sequence = "MIDS=" + df.sequence
    # print("MIDS")

df_sim = df[df.barcodeID.str.endswith("simulated")].reset_index(drop=True)
df_real = df[~df.barcodeID.str.endswith("simulated")].reset_index(drop=True)

# ====================================
# Define One-hot encording
# ====================================


class hot_dna:
    def __init__(self, fasta):
        if re.search(">", fasta):
            name = re.split("\n", fasta)[0]
            sequence = re.split("\n", fasta)[1]
        else:
            name = 'unknown_sequence'
            sequence = fasta
        seq_array = array(list(sequence))
        label_encoder = LabelEncoder()
        integer_encoded_seq = label_encoder.fit_transform(seq_array)
        # one hot the sequence
        onehot_encoder = OneHotEncoder(sparse=False, categories='auto')
        integer_encoded_seq = integer_encoded_seq.reshape(
            len(integer_encoded_seq), 1)
        onehot_encoded_seq = onehot_encoder.fit_transform(integer_encoded_seq)
        self.name = name
        self.sequence = fasta
        self.integer = integer_encoded_seq
        self.onehot = onehot_encoded_seq


def X_onehot(X_data):
    import collections
    char_num = np.zeros(10)
    for i in range(10):
        tmp = list(X_data.iloc[i][1])
        tmp = collections.Counter(tmp).values()
        char_num[i] = len(tmp)
    char_num = int(char_num.max())
    X = np.empty((len(X_data[0]), len(
        X_data.iloc[1, 1]), char_num), dtype="uint8")
    for i in tqdm(range(0, len(X_data[0]))):
        X[i] = hot_dna(X_data.iloc[i, 1]).onehot[0:len(X_data.iloc[1, 1])]
    X = X[:, :, 1:]
    return(X)

# ====================================
# Save One-hot matrix
# ====================================


print("One-hot encording simulated reads...")
X_sim = X_onehot(df_sim)
print("One-hot encording real reads...")
X_real = X_onehot(df_real)
np.savez_compressed(output_npz,
                    X_sim=X_sim,
                    X_real=X_real
                    )

# ====================================
# # Load One-hot matrix
# ====================================
np.load = partial(np.load, allow_pickle=True)
npz = np.load(output_npz)

X_sim = npz["X_sim"]
X_real = npz["X_real"]

X = X_sim
labels, labels_index = pd.factorize(df_sim.barcodeID)
labels_categorical = np_utils.to_categorical(labels)

X_train, X_test, Y_train, Y_test = train_test_split(
    X, labels_categorical,
    test_size=0.2, shuffle=True)

# ====================================
# # Model construction
# ====================================

tf.get_logger().setLevel('INFO')
model = Sequential()

model.add(Conv1D(filters=32, kernel_size=32, activation="relu",
                 input_shape=(X.shape[1], X.shape[2]), name="1st_Conv1D"))
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

#!#! L2 or FC
# L2 <<<<<<<
alpha = 0.1  # hyperparameter
model.add(Dense(64, activation='linear',
                activity_regularizer=regularizers.l2(alpha), name="L2-softmax"))

model.add(Dense(len(labels_index), activation='softmax', name="final_layer"))
model.compile(optimizer='adam', loss='categorical_crossentropy',
              metrics=['accuracy'])
model.summary()
# -

stack = model.fit(X_train, Y_train, epochs=20, verbose=0,
                  validation_split=0.2, shuffle=True)

# ====================================
# ## Save the model
# ====================================

model.save(output_model + '.h5')

# ====================================
# ## Plot loss and Accuracy
# ====================================
plt.figure()
plt.style.use('default')
plt.rcParams["font.size"] = 15
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['font.family'] = 'Arial'

fig = plt.figure(figsize=(15, 5))
ax1 = fig.add_subplot(1, 2, 1)
ax2 = fig.add_subplot(1, 2, 2)

ax1.plot(stack.history['loss'], lw=3)
ax1.plot(stack.history['val_loss'], lw=3)
ax1.set_title('Loss', fontsize=20)
# ax1.set_ylabel('loss')
ax1.set_xlabel('Epoch')
ax1.legend(['train', 'validation'])

ax2.plot(stack.history['accuracy'], lw=5)
ax2.plot(stack.history['val_accuracy'], lw=5)
ax2.set_title('Accuracy', fontsize=20)
# ax2.set_ylabel('accuracy')
ax2.set_xlabel('Epoch')
ax2.legend(['train', 'validation'])

# ------------------------
fig_name = "loss_acc"
for fig_dir in fig_dirs:
    fig_type = re.sub(".*/", "", fig_dir)
    plt.savefig(f"{fig_dir}/{output_figure}_{fig_name}.{fig_type}",
                dpi=350, bbox_inches="tight")


# ====================================
# ## Compute cosine similarity
# ====================================

def get_score_cosine(model, train, test):
    model = Model(model.get_layer(index=0).input,
                  model.get_layer(index=-2).output)  # Delete FC layer
    # print(model.summary())
    print("Obtain vectors from 1000 simulated reads...")
    normal_vector = model.predict(train, verbose=1, batch_size=1)
    print("Obtain vectors of each reads...")
    predict_vector = model.predict(test, verbose=1, batch_size=1)
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


X_all = np.concatenate((X_sim, X_real))
cos_all, normal_vector, predict_vector = get_score_cosine(
    model, X_train[0:1000], X_all)

df_name = pd.concat([df_sim[["barcodeID", "seqID"]].reset_index(drop=True),
                     df_real[["barcodeID", "seqID"]].reset_index(drop=True)])

df_all = pd.concat([df_name.reset_index(
    drop=True), pd.DataFrame(cos_all)], axis=1)
df_all.columns = ["barcodeID", "seqID", "cos_similarity"]

df_all["label"] = df_all.barcodeID.apply(
    lambda x: 1 if x.endswith("simulated") else -1)

optimal_threshold = df_all[df_all.label == 1].cos_similarity.quantile(0.001)

df_all["abnormal_prediction"] = df_all.cos_similarity.apply(
    lambda x: 1 if x > optimal_threshold else -1)

# df_unexpected = df_all[df_all.label == -1].reset_index()
# df_unexpected["abnormal_prediction"] = df_unexpected.cos_similarity.apply(
#     lambda x: "normal" if x > optimal_threshold else "abnormal")

# ====================================
# ## Plot cosine similarity
# ====================================

plt.figure(figsize=(6, 10))
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = '20'
sns.set_style("whitegrid")

ax = sns.boxplot(x="cos_similarity", y="barcode", data=df_all,
                 showfliers=False)
ax.set(xlim=[0.0, 1.05])
ax.set_title("Cosine similarity")
ax.set_xlabel("")
ax.set_ylabel("")
ax.axvline(x=optimal_threshold, ls="--", color="r",
           alpha=0.5, label="threshold of abnormality")
plt.legend(loc='upper left')
# ----------------------------
fig_name = "cosine_similarity"
for fig_dir in fig_dirs:
    fig_type = re.sub(".*/", "", fig_dir)
    plt.savefig(f"{fig_dir}/{output_figure}_{fig_name}.{fig_type}",
                dpi=350, bbox_inches="tight")

# ====================================
# ## Save anomaly classification for
# ## the detection of non-probrematic anomaly
# ====================================

output = df_real[["barcodeID", "seqID"]]
output["abnormal_prediction"] = df_all[df_all.label == -1].reset_index().cos_similarity.apply(
    lambda x: "Normal" if x > optimal_threshold else "Abnormal(problematic)")

output.to_csv(
    '.tmp_/anomaly_classification.txt',
    header=False, index=False, sep="\t")

pd.Series(labels_index).to_csv(
    '.tmp_/anomaly_classification_labels.txt',
    header=False, index=False, sep="\t")
