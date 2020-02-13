import numpy as np
import pandas as pd

import hdbscan

# ====================================
# # Import data
# ====================================
df = pd.read_csv(".tmp_/tmp_anomaly_classification",
                 header=None,
                 names=["seqID", "start", "cutsite", "cutlength"],
                 sep=' ')

X = df[["start", "cutsite", "cutlength"]]
X_array = np.array(X)

label = hdbscan.HDBSCAN(
    # min_samples=int(X_array.shape[0]/5),
    min_samples=10,
    min_cluster_size=int(X_array.shape[0]/5)+1,
).fit_predict(X_array)
pd.Series(label+1).value_counts()

df["label"] = label+2
df.to_csv(".tmp_/tmp_anomaly_label",
          sep="\t", header=False, index=False)
