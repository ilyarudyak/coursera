# Network Parameters
NUM_HIDDEN_UNITS = 256

# Training Parameters
BATCHES_PER_LOG = 100
CHECKPOINT_FILENAME = "checkpoints/default.checkpt"
LEARNING_RATE = 0.01
MEMORY_LEN = 1000000
MINIBATCH_SIZE = 32
NUM_BATCHES_PER_EPOCH = 50000
NUM_EPOCHS = 100

# Evaluation Parameters
NUM_EPOCHS_PER_EVAL = 1
NUM_POINTS_PER_EVALUATION = 1000
POINTS_PER_LOG = 100

# Game Params
NUM_ACTIONS = 2 # Can move up or down.
NUM_STATE_DIMENSIONS = 6 # 2x positions of paddle; pos and vel of ball

# Exploration Params
EXPLORATION_START_RATE = 0.90
EXPLORATION_DECAY_RATE = 0.01
EXPLORATION_RATE_MIN = 0.01
REWARD_DECAY_START = 0.0
REWARD_DECAY_GROWTH_FACTOR = 0.00
REWARD_DECAY_MAX = 0.0

TRAINING_MODE = False
CHOOSE_BEST_ALWAYS = False
TRAINING_MODE_NO_BOUNCES = False
REWARD_TYPE = "IDEAL_ANTICIPATION_REWARD"
REWARD_PROBABILITY = 1.0
REWARD_SCALING_FACTOR = 100
SCALE_REWARD_BY_DISTANCE_TO_PADDLE = True
DISTANCE_ERROR_POWER = 2.0
