#!/usr/bin/env python

from keras.models import Model
from keras.models import load_model as keras_load_model
from keras.layers import Dense, Dropout, Activation, Flatten, Conv1D, MaxPooling1D, Input, Lambda
from keras import metrics
from keras.callbacks import ModelCheckpoint, EarlyStopping
import tensorflow as tf
import argparse
import os, sys


parser = argparse.ArgumentParser(description = "Train DeepSea...)", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--train", nargs='+', help = "training data *.tfr generated by 'generate_tfr.py' ", required = True)
parser.add_argument("--val", nargs='+', help = "validation data *.tfr generated by 'generate_tfr.py", required = True)
parser.add_argument("--out", default = "deepsea_model", help = "model output files - only best model saved")
parser.add_argument('--batch_size', type=int, default = 1024, help = 'batch size for training')
parser.add_argument('--threads', type=int, default = 32, help = 'CPU cores for data pipeline loading')

args = parser.parse_args()

train_files = args.train
val_files = args.val
out = args.out
batch_size = args.batch_size
num_threads = args.threads

# train_deepsea.py --train *.tfr --val *.val --out deepsea_model.h5 --batch_size 1024


# singularity shell --nv /mnt/users/ngda/sofware/singularity/ndatth-deepsea-v0.0.0.img

# mv 21.tfr 21.val
# mv 25.tfr 25.val

# import glob
# batch_size = 1024
# num_threads = 32
# train_files = glob.glob('*.tfr')
# val_files = glob.glob('*.val')
# out = "deepsea_model"


out_model = out + '.h5'
out_hist = out + '.csv'

# Decoding function
def parse_record(record):
    name_to_features = {
        'seq': tf.io.FixedLenFeature([], tf.string),
        'label': tf.io.FixedLenFeature([], tf.string),
    }
    return tf.io.parse_single_example(record, name_to_features)



def decode_record(record):
    seq = tf.io.decode_raw(
        record['seq'], out_type=tf.float16, little_endian=True, fixed_length=None, name=None
    )
    label = tf.io.decode_raw(
        record['label'], out_type=tf.int8, little_endian=True, fixed_length=None, name=None
    )
    seq = tf.reshape(seq, [-1,4])
    #label = tf.cast(label, tf.float16)
    return (seq, label)



def get_dataset(record_file, num_threads = 8, batch = 512):
    dataset = tf.data.TFRecordDataset(record_file, num_parallel_reads = num_threads, compression_type = 'GZIP')
    dataset = dataset.map(parse_record, num_parallel_calls = num_threads)
    dataset = dataset.map(decode_record, num_parallel_calls = num_threads)
    dataset = dataset.shuffle(buffer_size = batch*10).batch(batch)
    return dataset
    

## source model from https://github.com/zj-zhang/deepsea_keras/blob/master/deepsea_keras/model.py

def build_model():
    inp = Input(shape=(NUM_INPUT, 4))
    x = Conv1D(filters=320, kernel_size=8, activation="relu", name="conv1d_1")(inp)
    x = MaxPooling1D(pool_size=4, strides=4)(x)
    x = Dropout(0.2)(x)
    x = Conv1D(filters=480, kernel_size=8, activation="relu", name="conv1d_2")(x)
    x = MaxPooling1D(pool_size=4, strides=4)(x)
    x = Dropout(0.2)(x)
    x = Conv1D(filters=960, kernel_size=8, activation="relu", name="conv1d_3")(x)
    x = Dropout(0.5)(x)
    x = Flatten(data_format="channels_first")(x)
    x = Dense(DENSE_UNITS, activation="relu", name="dense_1")(x)
    out = Dense(NUM_OUTPUT, activation="sigmoid", name="dense_2")(x)
    model = Model(inputs=[inp], outputs=[out])
    model.compile(
        optimizer = tf.keras.optimizers.Adam(learning_rate = 0.001),
        loss = tf.keras.losses.BinaryCrossentropy(),
        metrics = ['binary_accuracy', tf.keras.metrics.AUC(), tf.keras.metrics.Precision(), tf.keras.metrics.Recall()]
    )
    return model





## getting dims of input/output data:

data_dim = get_dataset(val_files, batch = 1)
dim_it = iter(data_dim)
dim_sample = next(dim_it)

DENSE_UNITS = 925
NUM_INPUT = dim_sample[0].shape[1]
NUM_OUTPUT = dim_sample[1].shape[1]


## model setting in multiple GPUs

# # Create a MirroredStrategy.
# strategy = tf.distribute.MirroredStrategy()
# print("Number of devices: {}".format(strategy.num_replicas_in_sync))

# with strategy.scope():
#     # Everything that creates variables should be under the strategy scope.
#     # In general this is only model construction & `compile()`.
#     model = build_model()


model = build_model()

checkpointer = ModelCheckpoint(filepath = out_model, verbose=1, save_best_only=True)
earlystopper = EarlyStopping(monitor="val_loss", patience=5, verbose=1)

train = get_dataset(train_files, batch= batch_size, num_threads = num_threads)
val = get_dataset(val_files, batch= batch_size, num_threads = num_threads)


history = model.fit(train, epochs=200, validation_data=val, callbacks=[checkpointer, earlystopper])


## export history
import pandas as pd
hist_df = pd.DataFrame(history.history) 
hist_df.to_csv(out_hist)


# # Plot utility
# def plot_graphs(hist_df, string):


#   plt.plot(history[string])
#   plt.plot(history['val_'+string])
#   plt.xlabel("Epochs")
#   plt.ylabel(string)
#   plt.legend([string, 'val_'+string])
#   plt.show()

# plot_graphs(hist_df, "accuracy")
# plot_graphs(hist_df, "loss")