import numpy as np
import os

BATCH_SIZES = [64, 128, 256, 512, 1024, 2048]
BATCHES_PER_LOGGING = [100, 50, 25, 10, 5, 2]
INITIAL_LEARNING_RATE = 0.001
KEEP_PROB = 0.50
LEARNING_RATE_DECAY = 0.80
MAX_ROTATION = 15 * ((2 * np.pi) / 360)
MAX_SHIFT = 2
MAX_ZOOM = 0.10
NUM_EPOCHS = len(BATCH_SIZES)
PROJECT_DIR = os.path.dirname(os.path.realpath(__file__))

DATA_DIR = f"{PROJECT_DIR}/data"