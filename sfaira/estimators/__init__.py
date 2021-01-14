from sfaira.estimators.keras import EstimatorKeras, EstimatorKerasEmbedding, EstimatorKerasCelltype

try:
    from sfaira_extension.estimators import *  # noqa: F403
except ImportError:
    pass
