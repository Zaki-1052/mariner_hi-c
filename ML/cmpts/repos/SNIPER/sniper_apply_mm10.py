# ML/cmpts/repos/SNIPER/sniper_apply_mm10.py
# mm10 application entry point for SNIPER.
# Inlines apply logic; imports mm10 utilities.
import os
import sys
import numpy as np

import tensorflow as tf
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

from scipy.io import loadmat, savemat
from keras.models import load_model

from ML.cmpts.repos.SNIPER.utilities.data_processing_mm10 import (
    hicToMat, trimMat, contactProbabilities, Sigmoid, predictionsToBed
)
from ML.cmpts.repos.SNIPER.utilities.input import get_application_params
from ML.cmpts.repos.SNIPER.utilities.mm10_config import (
    JUICER_CONTAINER, JUICER_JAR
)

JUICER_CMD = 'singularity exec --bind /scratch,/expanse {0} java -jar {1}'.format(
    JUICER_CONTAINER, JUICER_JAR
)


def apply_on_hic_mm10(params):
    inputM = hicToMat(params['input_file'], JUICER_CMD,
        tmp_dir=params['dump_dir'],
        prefix='input',
        autoremove=params['autoremove'],
        overwrite=params['overwrite'],
        save_matrix=params['save_matrix']
    )

    print('Trimming sparse, NA, and B4 regions...')
    inputM = trimMat(inputM, params['cropIndices'])
    print('Computing contact probabilities')
    inputM = contactProbabilities(inputM)

    odd_encoder = load_model(params['odd_encoder'])
    odd_clf = load_model(params['odd_classifier'])
    even_encoder = load_model(params['even_encoder'])
    even_clf = load_model(params['even_classifier'])

    odd_enc = Sigmoid(odd_encoder.predict(inputM))
    odd_predictions = odd_clf.predict(odd_enc)
    even_enc = Sigmoid(even_encoder.predict(inputM.T))
    even_predictions = even_clf.predict(even_enc)

    savemat(params['output_path'], {
        'odd_predictions'  : odd_predictions,
        'even_predictions' : even_predictions,
    })

    predictionsToBed(
        os.path.splitext(params['output_path'])[0] + '.bed',
        odd_predictions, even_predictions, params['cropMap']
    )


def apply_on_mat_mm10(params):
    inputM = loadmat(params['input_file'])['inter_matrix']

    print('Trimming sparse, NA, and B4 regions...')
    inputM = trimMat(inputM, params['cropIndices'])
    print('Computing contact probabilities')
    inputM = contactProbabilities(inputM)

    odd_encoder = load_model(params['odd_encoder'])
    odd_clf = load_model(params['odd_classifier'])
    even_encoder = load_model(params['even_encoder'])
    even_clf = load_model(params['even_classifier'])

    odd_enc = Sigmoid(odd_encoder.predict(inputM))
    odd_predictions = odd_clf.predict(odd_enc)
    even_enc = Sigmoid(even_encoder.predict(inputM.T))
    even_predictions = even_clf.predict(even_enc)

    savemat(params['output_path'], {
        'odd_predictions'  : odd_predictions,
        'even_predictions' : even_predictions,
    })


if __name__ == '__main__':
    params = get_application_params()

    params['cropMap'] = loadmat('crop_map/mm10_cropMap.mat')
    params['cropIndices'] = loadmat('crop_map/mm10_cropIndices.mat')

    if not params['usemat']:
        apply_on_hic_mm10(params)
    else:
        apply_on_mat_mm10(params)
