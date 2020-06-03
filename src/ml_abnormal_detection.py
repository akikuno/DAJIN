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
ont_cont = args[2]
mutation_type = args[3]
threads = int(args[4])

if ont_cont == "" :
    raise ValueError("ont_cont is empty")

if mutation_type == "" :
    raise ValueError("mutation_type is empty")

if threads == "" :
    threads = 1

# file_name = ".DAJIN_temp/data/DAJIN_MIDS.txt"
# ont_cont = "barcode01"
# threads = 10

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

# =============================================================
# Novelity detection
# =============================================================

def label_encorde_seq(seq):
    label_seq = seq.apply(list).\
        apply(LabelEncoder().fit_transform).\
        apply(pd.Series).\
        to_numpy()
    return(label_seq)


# =============================================================
# Train test split
# =============================================================

labels, labels_index = pd.factorize(df_sim.barcodeID)

# X_sim = label_encorde_seq(df_sim.seq)

X_train, X_test, Y_train, Y_test = train_test_split(
    label_encorde_seq(df_sim.seq), labels,
    test_size=0.2, shuffle=True)

del df_sim

# =============================================================
# Novelity detection
# =============================================================


clf = LocalOutlierFactor(n_neighbors=20, leaf_size = 400, novelty=True, n_jobs = threads)

clf.fit(X_train[:10000])

X_real = label_encorde_seq(df_real.seq)
del df_real["seq"]

outliers = clf.predict(X_real)
outliers = np.where(outliers==1, "normal", "abnormal") 
df_real["outliers"] = outliers

# df_real.groupby("barcodeID").outliers.value_counts()
df_real.to_csv(
    ".DAJIN_temp/data/DAJIN_MIDS_prediction_result.txt",
    header=False, index=False, sep="\t")

if mutation_type == "P" :
    sys.exit()

# =============================================================
# Classification
# =============================================================

model = lgb.LGBMClassifier(n_jobs=threads)
model.fit(X_train, Y_train)

prediction = model.predict_proba(X_real)
prediction = np.argmax(prediction, axis=1)

df_real["prediction"] = prediction
df_real["prediction"].mask(df_real["outliers"] == "abnormal", "abnormal", inplace=True)
del df_real["outliers"]

# df_real["prediction"].mask(df_real["outliers"] == 1, "normal", inplace=True)

for index,label in enumerate(labels_index):
    label=label.replace("_simulated","")
    df_real["prediction"].mask(df_real["prediction"] == index, label, inplace=True)


df_real.to_csv(
    ".DAJIN_temp/data/DAJIN_MIDS_prediction_result.txt",
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


# # Save labels
# pd.Series(labels_index).to_csv(
#     '.DAJIN_temp/data/DAJIN_anomaly_classification_labels.txt',
#     header=False, index=False, sep="\t")
