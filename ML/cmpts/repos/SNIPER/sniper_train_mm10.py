# ML/cmpts/repos/SNIPER/sniper_train_mm10.py
# mm10 training entry point for SNIPER.
# Inlines 4-class Classifier and training loop; imports mm10 utilities.
import os
import sys
import numpy as np

import tensorflow as tf
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

from scipy.io import loadmat
from keras.layers import Dense, Dropout
from keras.models import Sequential
from keras.utils import to_categorical

from ML.cmpts.repos.SNIPER.utilities.data_processing_mm10 import (
    hicToMat, trimMat, contactProbabilities, bootstrap, Sigmoid
)
from ML.cmpts.repos.SNIPER.pipeline.models import DenoisingAutoencoder
from ML.cmpts.repos.SNIPER.utilities.input import get_params
from ML.cmpts.repos.SNIPER.utilities.mm10_config import (
    N_CLASSES, TRAIN_VAL_SPLIT,
    JUICER_CONTAINER, JUICER_JAR
)

JUICER_CMD = 'singularity exec --bind /scratch,/expanse {0} java -jar {1}'.format(
    JUICER_CONTAINER, JUICER_JAR
)


def Classifier_mm10(X):
    input_dim = X.shape[1]

    clf = Sequential([
        Dense(64, activation='relu', input_dim=input_dim),
        Dropout(0.25),
        Dense(16, activation='relu'),
        Dropout(0.25),
        Dense(N_CLASSES, activation='softmax')
    ])

    clf.compile(loss='categorical_crossentropy', optimizer='rmsprop', metrics=['accuracy'])

    return clf


def trainNN_mm10(inputM, targetM, params):
    S = TRAIN_VAL_SPLIT
    assert inputM.shape[0] > S, (
        'Odd-chromosome bins ({0}) must exceed train/val split ({1}). '
        'Check crop map.'.format(inputM.shape[0], S)
    )
    assert inputM.shape[1] > S, (
        'Even-chromosome bins ({0}) must exceed train/val split ({1}). '
        'Check crop map.'.format(inputM.shape[1], S)
    )

    odd_labels = loadmat(params['label_file'])['rows'].flatten()
    even_labels = loadmat(params['label_file'])['cols'].flatten()

    print('Training autoencoders...')

    odd_dae_model, odd_encoder, _ = DenoisingAutoencoder(inputM, targetM)
    even_dae_model, even_encoder, _ = DenoisingAutoencoder(inputM.T, targetM.T)

    odd_dae_model.fit(inputM[:S], targetM[:S], epochs=10, batch_size=32,
            validation_data=[inputM[S:], targetM[S:]])
    even_dae_model.fit(inputM.T[:S], targetM.T[:S], epochs=10, batch_size=32,
            validation_data=[inputM.T[S:], targetM.T[S:]])

    odd_encodings = Sigmoid(odd_encoder.predict(inputM))
    even_encodings = Sigmoid(even_encoder.predict(inputM.T))

    odd_clf = Classifier_mm10(odd_encodings)
    even_clf = Classifier_mm10(even_encodings)

    odd_X, odd_y = bootstrap(odd_encodings[:S], odd_labels[:S], samplesPerClass=2000)
    even_X, even_y = bootstrap(even_encodings[:S], even_labels[:S], samplesPerClass=2000)

    print('Training classifiers...')
    odd_clf.fit(odd_X, to_categorical(odd_y, num_classes=N_CLASSES), epochs=25, batch_size=32, verbose=0,
        validation_data=[odd_encodings[S:], to_categorical(odd_labels[S:], num_classes=N_CLASSES)])

    even_clf.fit(even_X, to_categorical(even_y, num_classes=N_CLASSES), epochs=25, batch_size=32, verbose=0,
        validation_data=[even_encodings[S:], to_categorical(even_labels[S:], num_classes=N_CLASSES)])

    odd_dae_model.save(os.path.join(params['dump_dir'], 'odd_chrm_autoencoder.h5'))
    odd_encoder.save(os.path.join(params['dump_dir'], 'odd_chrm_encoder.h5'))
    odd_clf.save(os.path.join(params['dump_dir'], 'odd_chrm_classifier.h5'))

    even_dae_model.save(os.path.join(params['dump_dir'], 'even_chrm_autoencoder.h5'))
    even_encoder.save(os.path.join(params['dump_dir'], 'even_chrm_encoder.h5'))
    even_clf.save(os.path.join(params['dump_dir'], 'even_chrm_classifier.h5'))


def train_with_hic_mm10(params):
    print('Constructing input matrix')

    inputM = hicToMat(params['input_file'],
        JUICER_CMD,
        tmp_dir=params['dump_dir'],
        prefix='input',
        autoremove=params['autoremove'],
        overwrite=params['overwrite'],
        save_matrix=params['save_matrix']
    )

    print('Trimming sparse regions...')
    inputM = trimMat(inputM, params['cropIndices'])
    print('Computing contact probabilities')
    inputM = contactProbabilities(inputM)

    print('Constructing target matrix')

    targetM = hicToMat(params['target_file'],
        JUICER_CMD,
        tmp_dir=params['dump_dir'],
        prefix='target',
        autoremove=params['autoremove'],
        overwrite=params['overwrite'],
        save_matrix=params['save_matrix']
    )

    print('Trimming sparse regions...')
    targetM = trimMat(targetM, params['cropIndices'])
    print('Computing contact probabilities')
    targetM = contactProbabilities(targetM)

    trainNN_mm10(inputM, targetM, params)


def train_with_mat_mm10(params):
    print('Using pre-computed .mat files, skipping hic-to-mat')
    inputM = loadmat(params['input_file'])['inter_matrix']
    targetM = loadmat(params['target_file'])['inter_matrix']

    print('Trimming sparse regions from input matrix...')
    inputM = trimMat(inputM, params['cropIndices'])
    print('Computing contact probabilities')
    inputM = contactProbabilities(inputM)

    print('Trimming sparse regions from target matrix...')
    targetM = trimMat(targetM, params['cropIndices'])
    print('Computing contact probabilities')
    targetM = contactProbabilities(targetM)

    trainNN_mm10(inputM, targetM, params)


if __name__ == '__main__':
    params = get_params()

    params['cropMap'] = loadmat('crop_map/mm10_cropMap.mat')
    params['cropIndices'] = loadmat('crop_map/mm10_cropIndices.mat')

    if not params['usemat']:
        train_with_hic_mm10(params)
    else:
        train_with_mat_mm10(params)
