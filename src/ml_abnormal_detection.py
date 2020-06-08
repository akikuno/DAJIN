import os
os.environ['TF_CPP_MIN_LOG_LEVEL']='3'
os.environ['TF_FORCE_GPU_ALLOW_GROWTH']='true'

import sys
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder, LabelBinarizer
from sklearn.neighbors import LocalOutlierFactor

import tensorflow as tf
from tensorflow.keras.layers import Conv1D, Dense, Flatten,MaxPooling1D

################################################################################
#! I/O naming
################################################################################

#===========================================================
#? TEST auguments
#===========================================================

file_name = ".DAJIN_temp/data/DAJIN_MIDS.txt"
ont_cont = "barcode01"
mutation_type = "I"
threads = 10

#===========================================================
#? Auguments
#===========================================================

args = sys.argv
file_name = args[1]
ont_cont = args[2]
mutation_type = args[3]
threads = int(args[4])

if ont_cont == "" :
    raise ValueError("control is empty")

if mutation_type == "" :
    raise ValueError("mutation_type is empty")

if threads == "" :
    threads = 1

#===========================================================
#? Input
#===========================================================

df = pd.read_csv(file_name, header=None, sep='\t')
df.columns = ["seqID", "seq", "barcodeID"]
df.seq = "MIDS=" + df.seq

df_cont = df[df.barcodeID.str.contains(ont_cont)].reset_index(drop=True)
df_cont.barcodeID = "wt_simulated"

df_sim = df[df.barcodeID.str.contains("simulated")].reset_index(drop=True)
df_sim = df_sim.append(df_cont).reset_index(drop=True)
df_real = df[~df.barcodeID.str.contains("simulated")].reset_index(drop=True)

del df
del df_cont

################################################################################
#! Novelity (Anomaly) detection
################################################################################

#===========================================================
#? Encording
#===========================================================

def label_encorde_seq(seq):
    label_seq = seq.apply(list).\
        apply(LabelEncoder().fit_transform).\
        apply(pd.Series).\
        to_numpy()
    return(label_seq)

#===========================================================
#? Train test split
#===========================================================

labels, labels_index = pd.factorize(df_sim.barcodeID)

X_train, X_test, Y_train, Y_test = train_test_split(
    label_encorde_seq(df_sim.seq), labels,
    test_size=0.2, shuffle=True)

#===========================================================
#? LOF
#===========================================================
print("LOF") #? >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
clf = LocalOutlierFactor(n_neighbors=20, leaf_size = 400, novelty=True, n_jobs = threads)

clf.fit(X_train[:10000])

outliers = clf.predict(label_encorde_seq(df_real.seq))
outliers = np.where(outliers==1, "normal", "abnormal") 
df_real["outliers"] = outliers

# df_real.groupby("barcodeID").outliers.value_counts()
df_real.to_csv(
    ".DAJIN_temp/data/DAJIN_MIDS_prediction_result.txt",
    header=False, index=False, sep="\t")

if mutation_type == "P" :
    sys.exit()
    
################################################################################
#! Classification
################################################################################

#===========================================================
#? Encording
#===========================================================

def onehot_encode_seq(seq):
    onehot_seq = np.apply_along_axis(LabelBinarizer().fit_transform, 1, label_encorde_seq(seq))
    return(onehot_seq)

#===========================================================
#? Train test split
#===========================================================

labels = tf.keras.utils.to_categorical(labels)

X_train, X_test, Y_train, Y_test = train_test_split(
    onehot_encode_seq(df_sim.seq), labels,
    test_size=0.2, shuffle=True)

#===========================================================
#? CNN
#===========================================================

model = tf.keras.Sequential()

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

model.add(Dense(64, activation='relu', name="1st_FC"))

model.add(Dense(len(labels_index), activation='softmax', name="final_layer"))
model.compile(loss=tf.keras.losses.categorical_crossentropy,
              optimizer=tf.keras.optimizers.Adam(),
              metrics=['accuracy'])

model.fit(X_train, Y_train,
          epochs=15, verbose=1, batch_size = 64,
          validation_data=(X_test, Y_test), shuffle=True)

#===========================================================
#? prediction
#===========================================================
from tqdm import tqdm

X_real = onehot_encode_seq(df_real.seq)
del df_real["seq"] #! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

iter_ = 1000
prediction = np.zeros(X_real.shape[0], dtype="uint8")
for i in tqdm(range(0, X_real.shape[0], iter_)):
    predict_ = model.predict(X_real[i: i + iter_].astype("float16"),
                             verbose=0, batch_size=32)
    prediction[i:i+iter_] = np.argmax(predict_, axis=1)

del X_real #! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
# prediction = model.predict(onehot_encode_seq(df_real.seq).astype(np.float16))
# prediction = prediction.argmax(axis=1)

df_real["prediction"] = prediction
df_real["prediction"].mask(df_real["outliers"] == "abnormal", "abnormal", inplace=True)
del df_real["outliers"] #! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

for index,label in enumerate(labels_index):
    label=label.replace("_simulated","")
    df_real["prediction"].mask(df_real["prediction"] == index, label, inplace=True)

print(df_real.groupby("barcodeID").prediction.value_counts()) #? >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

df_real.to_csv(
    ".DAJIN_temp/data/DAJIN_MIDS_prediction_result.txt",
    header=False, index=False, sep="\t")




# # ====================================
# # RandomForestClassifier
# # ====================================
# from sklearn.ensemble import RandomForestClassifier
# clf = RandomForestClassifier(n_jobs = threads)
# clf.fit(X_train, Y_train)

# # test = df_real.seqID == "00222820-ad67-4774-b5f9-8fda0fba4c83"
# # test = df_real.seqID == "00ea11a6-9c44-4d42-84e7-e03a113c6cea"
# # test = X_real[test]
# # clf.predict(test)

# prediction = clf.predict(X_real)

# df_real["outliers"] = outliers
# df_real["prediction"] = prediction
# df_real["prediction"].mask(df_real["outliers"] == "abnormal", "abnormal", inplace=True)
# del df_real["outliers"] #!--------------------------------

# for index,label in enumerate(labels_index):
#     label=label.replace("_simulated","")
#     df_real["prediction"].mask(df_real["prediction"] == index, label, inplace=True)

# df_real.groupby("barcodeID").prediction.value_counts()

# df_real.to_csv(
#     ".DAJIN_temp/data/DAJIN_MIDS_prediction_result.txt",
#     header=False, index=False, sep="\t")


# =============================================================
# PCA
# =============================================================
# X_all = label_encorde_seq(df.seq)
# X_all = onehot_seq

# import numpy as np
# from sklearn.decomposition import PCA

# pca = PCA(n_components=5)
# X_all = pca.fit(X_all).transform(X_all)
# X_all = X_all * pca.explained_variance_ratio_

# X_sim = X_all[df.barcodeID.str.contains(ont_cont) + df.barcodeID.str.contains("simulated")]
# X_real = X_all[~df.barcodeID.str.contains("simulated")]

# nsamples, nx, ny = X_all.shape
# X_all_reshape = X_all.reshape((nsamples,nx*ny))
# X_sim = X_all_reshape[df.barcodeID.str.contains(ont_cont) | df.barcodeID.str.contains("simulated")]
# X_real = X_all_reshape[~df.barcodeID.str.contains("simulated")]


# =============================================================
# LGBMClassifier
# =============================================================
# import lightgbm as lgb
# from sklearn.metrics import accuracy_score

# model = lgb.LGBMClassifier(n_jobs=threads)
# model.fit(X_train, Y_train)

# prediction = model.predict_proba(X_test)
# prediction = np.argmax(prediction, axis=1)

# accuracy = accuracy_score(prediction, Y_test)
# print('Averaging accuracy:', accuracy)

# prediction = model.predict_proba(X_real)
# prediction = np.argmax(prediction, axis=1)

# df_real["prediction"] = prediction
# df_real["prediction"].mask(df_real["outliers"] == "abnormal", "abnormal", inplace=True)
# del df_real["outliers"]

# for index,label in enumerate(labels_index):
#     label=label.replace("_simulated","")
#     df_real["prediction"].mask(df_real["prediction"] == index, label, inplace=True)


# df_real.groupby("barcodeID").prediction.value_counts()

# df_real.to_csv(
#     ".DAJIN_temp/data/DAJIN_MIDS_prediction_result.txt",
#     header=False, index=False, sep="\t")

# ====================================
# IsolationForest
# ====================================
# from sklearn.ensemble import IsolationForest
# clf = IsolationForest(n_jobs=10) 

# train = X_train[0:10000]
# clf.fit(train)

# y_pred_outliers = clf.predict(X_real)

# df_real["outliers"] = y_pred_outliers

# df_real.groupby("barcodeID").outliers.value_counts()


# # Save labels
# pd.Series(labels_index).to_csv(
#     '.DAJIN_temp/data/DAJIN_anomaly_classification_labels.txt',
#     header=False, index=False, sep="\t")



# ====================================
# PCA for X_real
# ====================================
# import numpy as np
# from sklearn.decomposition import PCA

# pca = PCA(n_components=10)
# X_real = pca.fit(X_real).transform(X_real)
# print(pca.explained_variance_ratio_)

# colors = ['navy', 'turquoise', 'darkorange']
# lw = 2

# y, barcode = pd.factorize(df_real.barcodeID)

# import matplotlib.pyplot as plt
# plt.figure()
# for color, i, bar in zip(colors, [0, 1, 2], barcode):
#     plt.scatter(X_r[y == i, 0], X_r[y == i, 1], color=color, alpha=.8, lw=lw, label=bar)

# plt.legend(loc='best', shadow=False, scatterpoints=1)
# plt.savefig("test.png")

# ====================================
# Neural net
# ====================================

# from sklearn.neural_network import MLPClassifier
# clf = MLPClassifier(solver='lbfgs', hidden_layer_sizes=(5,2))
# clf.fit(X_train, Y_train)
# clf.score(X_test, Y_test)

# prediction = clf.predict_proba(X_real)
# prediction = np.argmax(prediction, axis=1)

# df_real["outliers"] = outliers
# df_real["prediction"] = prediction
# df_real["prediction"].mask(df_real["outliers"] == "abnormal", "abnormal", inplace=True)

# for index,label in enumerate(labels_index):
#     label=label.replace("_simulated","")
#     df_real["prediction"].mask(df_real["prediction"] == index, label, inplace=True)


# df_real.groupby("barcodeID").prediction.value_counts()
