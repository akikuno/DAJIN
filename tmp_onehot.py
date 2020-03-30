import numpy as np
from numpy import array
from numpy import argmax
from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import OneHotEncoder
import tensorflow as tf
# define example
# data = ['cold', 'cold', 'warm', 'cold', 'hot', 'hot', 'warm', 'cold', 'warm', 'hot']
# data = "ACTG"

def one_hot_encording(seq):
    values = [i for i in seq]
    # values = np.char.split(values, sep ='') 
    #print(values)
    # integer encode
    label_encoder = LabelEncoder()
    integer_encoded = label_encoder.fit_transform(values)
    # print(integer_encoded)
    # binary encode
    onehot_encoder = OneHotEncoder(sparse=False, dtype="uint8")
    integer_encoded = integer_encoded.reshape(len(integer_encoded), 1)
    onehot_encoded = onehot_encoder.fit_transform(integer_encoded)
    return onehot_encoded


def convert_sparse_matrix_to_sparse_tensor(X):
    coo = X.tocoo()
    indices = np.mat([coo.row, coo.col]).transpose()
    return tf.SparseTensor(indices, coo.data, coo.shape)

data = df.seq
X_all = np.zeros([len(data),len(data[0]),4], dtype="uint8")
# for i in len(data): 
for i in range(1):
    X_all[i,:,:] = one_hot_encording(data)


indice = convert_sparse_matrix_to_sparse_tensor(hoge)


# ====================================
# Define One-hot encording
# ====================================


class hot_dna:
    def __init__(self, fasta):
        seq_array = array(list(sequence))
        label_encoder = LabelEncoder()
        integer_encoded_seq = label_encoder.fit_transform(seq_array)
        # one hot the sequence
        onehot_encoder = OneHotEncoder(sparse=False, categories='auto', dtype='uint8')
        integer_encoded_seq = integer_encoded_seq.reshape(
            len(integer_encoded_seq), 1)
        onehot_encoded_seq = onehot_encoder.fit_transform(integer_encoded_seq)
        self.name = name
        self.sequence = fasta
        self.integer = integer_encoded_seq
        self.onehot = onehot_encoded_seq


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

# ====================================
# Save One-hot matrix
# ====================================


print("One-hot encording simulated reads...")
X_sim = X_onehot(df_sim.seq)
print("One-hot encording real reads...")
X_real = X_onehot(df_real.seq)
np.savez_compressed(output_npz,
                    X_sim=X_sim,
                    X_real=X_real
                    )
