from sfaira.estimators.keras.base import EstimatorKeras, EstimatorKerasEmbedding, EstimatorKerasCelltype

try:
    from sfaira_extension.estimators import *  # noqa: F403
except ImportError:
    pass
