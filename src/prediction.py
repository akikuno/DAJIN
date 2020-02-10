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
# # Import data
# ====================================

args = sys.argv
file_name = args[1]
df_anomaly = pd.read_csv(".tmp_/anomaly_classification_revised.txt",
                      header=None, sep='\t')

fig_dirs = ["results/figures/png", "results/figures/svg"]

output_figure = file_name.replace(".txt.gz", "").replace("data_for_ml/", "")
output_model = file_name.replace(".txt.gz", "").replace(
    "data_for_ml", "data_for_ml/model")

# ====================================
# # Load One-hot matrix
# ====================================
np.load = partial(np.load, allow_pickle=True)
npz = np.load(output_npz)

X_sim = npz["X_sim"]
X_real = npz["X_real"]

X = X_sim
labels, labels_id = pd.factorize(df_sim.iloc[:, 0])
labels_categorical = np_utils.to_categorical(labels)

X_train, X_test, Y_train, Y_test = train_test_split(
    X, labels_categorical,
    test_size=0.2, shuffle=True)

# ====================================
# # Load model
# ====================================

tf.get_logger().setLevel('INFO')
model = load_model(output_model + '.h5')

# ====================================
# # Prediction
# ====================================

print("Predict labels...")
predict = model.predict(X_real, verbose=1, batch_size=1)

predict = np.argmax(predict, axis=1)

df_predict = pd.Series(predict, dtype="str") + "_"

for i, j in enumerate(pd.Series(labels_id).str.replace("_simulated", "")):
    df_predict = df_predict.str.replace(str(i)+"_", j)

df_result = pd.DataFrame({"barcodeID": df_real.iloc[:, 0],
                          "seqID": df_real.iloc[:, 2],
                          "predict": df_predict,
                          "anomaly": df_unexpected.abnormal_prediction})

df_result.predict = df_result.predict.where(
    df_result.anomaly == "normal", "abnormal")
