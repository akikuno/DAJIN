# import warnings
# warnings.filterwarnings("ignore")
# warnings.simplefilter(action='ignore', category=FutureWarning)

# import re

# import matplotlib.pyplot as plt
# import seaborn as sns
import sys
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder

from sklearn.neighbors import LocalOutlierFactor
from sklearn.ensemble import IsolationForest

import lightgbm as lgb

# ====================================
# Input and format data
# ====================================

args = sys.argv
file_name = args[1]
control = args[2]
if control == "" :
    raise ValueError("control is empty")
# file_name = ".DAJIN_temp/data/DAJIN_MIDS.txt"
# control = "barcode01"
df = pd.read_csv(file_name, header=None, sep='\t')
df.columns = ["seqID", "seq", "barcodeID"]
df.seq = "MIDS=" + df.seq

df_sim = df[df.barcodeID.str.contains("simulated")].reset_index(drop=True)
df_cont = df[df.barcodeID.str.contains(control)].reset_index(drop=True)
df_cont.barcodeID = "wt_simulated"
df_sim = df_sim.append(df_cont).reset_index(drop=True)
df_real = df[~df.barcodeID.str.contains("simulated")].reset_index(drop=True)

del df
del df_cont

# ====================================
# Preporcessing...
# ====================================

def label_encorde_seq(seq):
    label_encoder = LabelEncoder()
    label_seq = seq.apply(lambda x: list(x)).\
        apply(lambda x: label_encoder.fit_transform(x)).\
        apply(lambda x: pd.Series(x)).\
        to_numpy()
    return(label_seq)

X_sim = label_encorde_seq(df_sim.seq)

# # ====================================
# # # Train test split...
# # ====================================

# print("Model training...")

labels, labels_index = pd.factorize(df_sim.barcodeID)
X_train, X_test, Y_train, Y_test = train_test_split(
    X_sim, labels,
    test_size=0.2, shuffle=True)


# ====================================
# LOF
# ====================================

clf = LocalOutlierFactor(n_neighbors=20, leaf_size = 400, novelty=True, n_jobs = 10)

train = X_train[:10000]
clf.fit(train)

del X_sim
X_real = label_encorde_seq(df_real.seq)

y_pred_outliers = clf.predict(X_real)
df_real["outliers"] = y_pred_outliers

df_real.groupby("barcodeID").outliers.value_counts()


# # ====================================
# # # lightGBM
# # ====================================

model = lgb.LGBMClassifier(n_jobs=10)
model.fit(X_train, Y_train)

y_pred = model.predict_proba(X_test)
y_pred_max = np.argmax(y_pred, axis=1)

accuracy = sum(Y_test == y_pred_max) / len(Y_test)
print(accuracy)

y_pred = model.predict_proba(X_real)
y_pred_max = np.argmax(y_pred, axis=1)

df_real["prediction"].mask(df_real["outliers"] == -1, "abnormal", inplace=True)
# df_real["prediction"].mask(df_real["outliers"] == 1, "normal", inplace=True)

for index,label in enumerate(labels_index):
    print(index)
    label=label.replace("_simulated","")
    df_real["prediction"].mask(df_real["prediction"] == index, label, inplace=True)

df_real.groupby("barcodeID").prediction.value_counts()



#! ====================================================

df_anomaly.to_csv(
    '.DAJIN_temp/data/DAJIN_anomaly_classification.txt',
    header=False, index=False, sep="\t")

# Save labels
pd.Series(labels_index).to_csv(
    '.DAJIN_temp/data/DAJIN_anomaly_classification_labels.txt',
    header=False, index=False, sep="\t")

# # # ====================================
# # # # IsolationForest
# # # ====================================
# clf = IsolationForest(n_jobs=10) 

# train = X_train[0:10000]
# clf.fit(train)

# y_pred_outliers = clf.predict(X_real)

# df_real["outliers"] = y_pred_outliers

# df_real.groupby("barcodeID").outliers.value_counts()

