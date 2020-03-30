import os
os.environ['TF_CPP_MIN_LOG_LEVEL']='3'
import warnings
warnings.filterwarnings("ignore")
warnings.simplefilter(action='ignore', category=FutureWarning)

import sys
import re
from functools import partial

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.model_selection import train_test_split
from tqdm import tqdm

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
df.seq = "MIDS=" + df.seq
df.columns = ["seqID", "seq", "barcodeID"]

df_sim = df[df.barcodeID.str.contains("simulated")].reset_index(drop=True)
df_real = df[~df.barcodeID.str.contains("simulated")].reset_index(drop=True)

# Output names
fig_dirs = ["results/figures/png", "results/figures/svg"]
output_npz = file_name.replace(".txt", ".npz").replace(
    "data_for_ml/", "data_for_ml/model/")
output_figure = file_name.replace(".txt", "").replace("data_for_ml/", "")
output_model = file_name.replace(".txt", "").replace(
    "data_for_ml", "data_for_ml/model")

# # ====================================
# # One-hot encording
# # ====================================
def X_onehot(X_data):
    X = np.empty((len(X_data), len(X_data[0]), 5),
        dtype="uint8")
    for i in tqdm(range(0, len(X_data))):
        sequence = X_data[i]
        seq_array = array(list(sequence))
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
    print("Obtain vectors from 1000 simulated reads...")
    normal_vector = model_.predict(train, verbose=1, batch_size=32)
    print("Obtain vectors of each reads...")
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
    lambda x: 1 if x > optimal_threshold else -1)
df_all = df_all._to_pandas()

# ====================================
# Output the results
# ====================================
# Save One-hot matrix
np.savez_compressed(output_npz,
                X_sim=X_sim,
                X_real=X_real
                )
# Save the model
model.save(output_model + '.h5')

# Save Anomaly annotation
df_anomaly = df_real[["barcodeID", "seqID"]]
df_anomaly["abnormal_prediction"] = df_all[df_all.label == -1].reset_index().cos_similarity.apply(
    lambda x: "normal" if x > optimal_threshold else "abnormal")
df_anomaly.to_csv(
    '.tmp_/DAJIN_anomaly_classification.txt',
    header=False, index=False, sep="\t")

# Save labels
pd.Series(labels_index).to_csv(
    '.tmp_/DAJIN_anomaly_classification_labels.txt',
    header=False, index=False, sep="\t")



# ====================================
# # ## Plot loss and Accuracy
# # ====================================
# plt.figure()
# plt.style.use('default')
# plt.rcParams["font.size"] = 15
# plt.rcParams['axes.linewidth'] = 1.5
# plt.rcParams['font.family'] = 'Arial'

# fig = plt.figure(figsize=(15, 5))
# ax1 = fig.add_subplot(1, 2, 1)
# ax2 = fig.add_subplot(1, 2, 2)

# ax1.plot(stack.history['loss'], lw=3)
# ax1.plot(stack.history['val_loss'], lw=3)
# ax1.set_title('Loss', fontsize=20)
# # ax1.set_ylabel('loss')
# ax1.set_xlabel('Epoch')
# ax1.legend(['train', 'validation'])

# ax2.plot(stack.history['accuracy'], lw=5)
# ax2.plot(stack.history['val_accuracy'], lw=5)
# ax2.set_title('Accuracy', fontsize=20)
# # ax2.set_ylabel('accuracy')
# ax2.set_xlabel('Epoch')
# ax2.legend(['train', 'validation'])

# # ------------------------
# fig_name = "loss_acc"
# for fig_dir in fig_dirs:
#     fig_type = re.sub(".*/", "", fig_dir)
#     plt.savefig(f"{fig_dir}/{output_figure}_{fig_name}.{fig_type}",
#                 dpi=350, bbox_inches="tight")


# # ====================================
# # ## Plot cosine similarity
# # ====================================

# plt.figure(figsize=(6, 10))
# plt.rcParams['axes.linewidth'] = 1.5
# plt.rcParams['font.family'] = 'Arial'
# plt.rcParams['font.size'] = '20'
# sns.set_style("whitegrid")

# ax = sns.boxplot(x="cos_similarity", y="barcodeID", data=df_all.sort_values(by="barcodeID"),
#                 showfliers=False)
# ax.set(xlim=[0.0, 1.05])
# ax.set_title("Abnormal allele detection")
# ax.set_xlabel("Cosine similarity")
# ax.set_ylabel("")
# ax.axvline(x=optimal_threshold, ls="--", color="r",
#         alpha=0.5, label="threshold of abnormality")
# plt.legend(loc='bottom left')
# # ----------------------------
# fig_name = "cosine_similarity"
# for fig_dir in fig_dirs:
#     fig_type = re.sub(".*/", "", fig_dir)
#     plt.savefig(f"{fig_dir}/{output_figure}_{fig_name}.{fig_type}",
#                 dpi=350, bbox_inches="tight")
