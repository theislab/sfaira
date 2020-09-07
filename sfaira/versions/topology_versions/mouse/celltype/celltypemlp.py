CELLTYPEMLP_TOPOLOGIES = {
    "0.0.1": {
        "genome": "Mus_musculus_GRCm38_97",
        "hyper_parameters": {
            "units": [],
            "activation": None,
            "use_bias": True,
            "l1_coef": 0.,
            "l2_coef": 0.,
            "kernel_initializer": 'glorot_uniform',
            "bias_initializer": 'zeros',
            "bias_regularizer": None,
            "activity_regularizer": None,
            "kernel_constraint": None,
            "bias_constraint": None
        }
    },
    "0.1.1": {
        "genome": "Mus_musculus_GRCm38_97",
        "hyper_parameters": {
            "units": [128],
            "activation": "relu",
            "use_bias": True,
            "l1_coef": 0.,
            "l2_coef": 0.,
            "kernel_initializer": 'glorot_uniform',
            "bias_initializer": 'zeros',
            "bias_regularizer": None,
            "activity_regularizer": None,
            "kernel_constraint": None,
            "bias_constraint": None
        }
    },
    "0.1.2": {
        "genome": "Mus_musculus_GRCm38_97",
        "hyper_parameters": {
            "units": [256, 128],
            "activation": "relu",
            "use_bias": True,
            "l1_coef": 0.,
            "l2_coef": 0.,
            "kernel_initializer": 'glorot_uniform',
            "bias_initializer": 'zeros',
            "bias_regularizer": None,
            "activity_regularizer": None,
            "kernel_constraint": None,
            "bias_constraint": None
        }
    },
    "0.1.3": {
        "genome": "Mus_musculus_GRCm38_97",
        "hyper_parameters": {
            "units": [512, 256, 128],
            "activation": "relu",
            "use_bias": True,
            "l1_coef": 0.,
            "l2_coef": 0.,
            "kernel_initializer": 'glorot_uniform',
            "bias_initializer": 'zeros',
            "bias_regularizer": None,
            "activity_regularizer": None,
            "kernel_constraint": None,
            "bias_constraint": None
        }
    }
}

# Load versions from extension if available:
try:
    import sfaira_extension as sfairae
    ADD_TOPOLOGIES = sfairae.versions.topology_versions.mouse.celltype.CELLTYPEMLP_TOPOLOGIES
    for k in CELLTYPEMLP_TOPOLOGIES.keys():
        if k in ADD_TOPOLOGIES.keys():
            CELLTYPEMLP_TOPOLOGIES.update(ADD_TOPOLOGIES)
except ImportError:
    pass
