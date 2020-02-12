import hdbscan
import time
from scipy.sparse.csgraph import connected_components
import umap
from sklearn.preprocessing import StandardScaler
from sklearn import cluster
from matplotlib.axes._axes import _log as matplotlib_axes_logger
import sys
import itertools
from sklearn.metrics import confusion_matrix
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.model_selection import train_test_split
from functools import partial
import os
import numpy as np
import pandas as pd
import re
from sklearn.preprocessing import OneHotEncoder
from sklearn.preprocessing import LabelEncoder
from sklearn.cluster import DBSCAN
from sklearn import metrics
from numpy import argmax
from numpy import array
import warnings
from tqdm import tqdm
warnings.filterwarnings('ignore')


# ====================================
# # Import data
# ====================================
testdata = pd.read_csv("tmp_testdata2",
                       header=None,
                       names=["seqID","start","cutsite","cutlength"],
                       sep=' ')

# file_name = "data_for_ml/sequence_MIDS.txt.gz"
# output_npz = file_name.replace(".txt.gz", ".npz").replace(
#     "data_for_ml/", "data_for_ml/model/")

# np.load = partial(np.load, allow_pickle=True)
# npz = np.load(output_npz)
# X_real = npz["X_real"]

# df_anomaly = pd.read_csv(".tmp_/anomaly_classification_revised.txt",
#                          header=None,
#                          names=["barcodeID", "seqID", "label"],
#                          sep='\t')

# labels_index = pd.read_csv(
#     ".tmp_/anomaly_classification_labels.txt",
#     header=None, sep='\t')

fig_dirs = ["results/figures/png", "results/figures/svg"]
output_figure = file_name.replace(".txt.gz", "").replace("data_for_ml/", "")
output_model = file_name.replace(".txt.gz", "").replace(
    "data_for_ml", "data_for_ml/model")

# ====================================
# # Extract barcode 12 and target_deletion
# ====================================
# testdata_list = tuple(testdata.seqID.to_list())
# objects = np.array(df_anomaly.barcodeID == "barcode12") * \
#     np.array(df_anomaly.label == "Abnormal(target_deletion)") * \
#         np.array(df_anomaly.seqID.str.startswith(testdata_list))

# df_b12 = df_anomaly[objects].reset_index(drop=True)
# df_b12["true_label"] = testdata.label

# X_b12 = X_real[objects]
# X = X_b12
X = testdata[["start","cutsite","cutlength"]]
X_reshape = np.array(X)
# -------------------------------------------------------
#df_b12 = df_b12[0:1000]
#X = X_b12[0:1000]
# -------------------------------------------------------
# X.shape
# X_reshape = X.reshape(X.shape[0], X.shape[1]*X.shape[2])
# X_reshape.shape

# np.seterr(divide='ignore', invalid='ignore')


# ====================================
# PCA/UMAP to HDBSCAN
# ====================================
# start = time.time()
# comp = 3
# X_pca = StandardScaler().fit_transform(X_reshape)
# pca = PCA(n_components=comp)
# X_pca = pca.fit_transform(X_pca)
# pca_time = time.time() - start
# print("pca_time:{0}".format(pca_time) + "[sec]")

# start = time.time()
# X_umap = umap.UMAP(n_neighbors=30,
#                    min_dist=0.0,
#                    n_components=2).fit_transform(X_reshape)
# umap_time = time.time() - start
# print("umap_time:{0}".format(umap_time) + "[sec]")

# db = cluster.DBSCAN(eps=3, min_samples=20, n_jobs=-1).fit_predict(X_umap)
# db = cluster.DBSCAN(eps=1000, min_samples=20, n_jobs=-1).fit_predict(X_pca)
# label = db
# label = testdata.label

inputdata = X_reshape
label = hdbscan.HDBSCAN(
    min_samples=int(inputdata.shape[0]/5),
    min_cluster_size=int(inputdata.shape[0]/5)+1,
).fit_predict(inputdata)
pd.Series(label).value_counts()

testdata["label"] = label
testdata.to_csv("hoge.txt", sep="\t", header=False, index=False)

umap_embedding = umap.UMAP().fit_transform(inputdata)
result_UMAP = pd.DataFrame(umap_embedding, columns=[
    'dim%i' % i for i in range(2)])
result_UMAP["label"] = label


fig = plt.figure(figsize=(8, 6), dpi=350)
ax = fig.add_subplot(1, 1, 1)

for i in np.unique(label):
    print(i)
    ax.scatter(result_UMAP[result_UMAP.label == i]["dim0"],
               result_UMAP[result_UMAP.label == i]["dim1"],
               # c=colorlist[i],
               c=plt.cm.tab20(i),
               alpha=0.5)

plt.title('UMAP of X_reshape')

fig_name = "umap_test"
for fig_dir in fig_dirs:
    fig_type = re.sub(".*/", "", fig_dir)
    plt.savefig(f"{fig_dir}/{output_figure}_{fig_name}.{fig_type}",
                dpi=350, bbox_inches="tight")


# op = cluster.OPTICS(min_samples=3, n_jobs=-1).fit_predict(X_pca)
# op
# km = cluster.KMeans(n_clusters=4,
#                     init='k-means++',
#                     max_iter=300,
#                     n_init=10,
#                     n_jobs=-1).fit_predict(X_pca)
# km

# ====================================================
# Run PCA
# ====================================================
label = df_b12.true_label
pca = PCA(n_components=3)
pca.fit(X_reshape)

# Store results of PCA in a data frame
result_PCA = pd.DataFrame(pca.transform(X_reshape), columns=[
                          'dim%i' % i for i in range(3)])
result_PCA["label"] = label
colorlist = ['#1f77b4', '#d62728', '#4daf4a', '#9467bd', '#999999']

# %matplotlib inline
plt.style.use('default')
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = '20'
matplotlib_axes_logger.setLevel('ERROR')

fig = plt.figure(figsize=(8, 6), dpi=350)
ax = fig.add_subplot(1, 1, 1)

for i in np.unique(label):
    print(i)
    ax.scatter(result_PCA[result_PCA.label == i]["dim0"],
               result_PCA[result_PCA.label == i]["dim1"],
               c=colorlist[i],
               alpha=0.5)

# plt.legend(bbox_to_anchor=(1.01, 1), loc='upper left',
#            borderaxespad=0, fontsize=10)
plt.title('PCA of X_pca')

fig_name = "pca_test"
for fig_dir in fig_dirs:
    fig_type = re.sub(".*/", "", fig_dir)
    plt.savefig(f"{fig_dir}/{output_figure}_{fig_name}.{fig_type}",
                dpi=350, bbox_inches="tight")


# # ====================================
# # Run UMAP
# # ====================================
# label = db

umap_embedding = umap.UMAP().fit_transform(X_reshape)
result_UMAP = pd.DataFrame(umap_embedding, columns=[
    'dim%i' % i for i in range(2)])
#
result_UMAP["label"] = label

fig = plt.figure(figsize=(8, 6), dpi=350)
ax = fig.add_subplot(1, 1, 1)

for i in np.unique(label):
    print(i)
    ax.scatter(result_UMAP[result_UMAP.label == i]["dim0"],
               result_UMAP[result_UMAP.label == i]["dim1"],
               #c=colorlist[i],
               c=plt.cm.tab20(i),
               alpha=0.5)

plt.title('UMAP of X_reshape')

fig_name = "umap_test"
for fig_dir in fig_dirs:
    fig_type = re.sub(".*/", "", fig_dir)
    plt.savefig(f"{fig_dir}/{output_figure}_{fig_name}.{fig_type}",
                dpi=350, bbox_inches="tight")


# # ==================================
# # # Elbow
# # ==================================

# distortions = []

# for i in range(1, 11):
#     km = cluster.KMeans(n_clusters=i,
#                         init='k-means++',
#                         n_init=10,
#                         max_iter=300,
#                         random_state=0,
#                         n_jobs=-1)
#     km.fit(X_reshape)
#     distortions.append(km.inertia_)

# plt.figure(figsize=(8, 6), dpi=350)
# plt.plot(range(1, 11), distortions, marker='o')
# plt.xlabel('Number of clusters')
# plt.ylabel('Distortion')
# plt.show()

# fig_name = "pca_test_elbow"
# for fig_dir in fig_dirs:
#     fig_type = re.sub(".*/", "", fig_dir)
#     plt.savefig(f"{fig_dir}/{output_figure}_{fig_name}.{fig_type}",
#                 dpi=350, bbox_inches="tight")


# # ====================================
# # Compute DBSCAN
# # ====================================
# #X_reshape = StandardScaler().fit_transform(X_reshape)

# db = cluster.DBSCAN(eps=10, min_samples=2, n_jobs=4).fit_predict(X_reshape)
# db
# op = cluster.OPTICS(min_samples=3, n_jobs=-1).fit_predict(X_reshape)
# op
# km = cluster.KMeans(n_clusters=4,
#                     init='k-means++',
#                     max_iter=300,
#                     n_init=10,
#                     n_jobs=-1).fit_predict(X_reshape)
# km
