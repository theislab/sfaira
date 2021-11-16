CELLTYPEMLP_TOPOLOGIES = {
    "0.0.1": {
        "model_type": "mlp",
        "input": {
            "genome": "Homo_sapiens.GRCh38.102",
            "genes": ["biotype", "protein_coding"],
        },
        "output": {
            "cl": "v2021-02-01",
            "targets": None,
        },
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
        "model_type": "mlp",
        "input": {
            "genome": "Homo_sapiens.GRCh38.102",
            "genes": ["biotype", "protein_coding"],
        },
        "output": {
            "cl": "v2021-02-01",
            "targets": None,
        },
        "hyper_parameters": {
            "units": [128],
            "activation": "selu",
            "use_bias": True,
            "l1_coef": 0.,
            "l2_coef": 0.,
            "kernel_initializer": 'lecun_normal',
            "bias_initializer": 'zeros',
            "bias_regularizer": None,
            "activity_regularizer": None,
            "kernel_constraint": None,
            "bias_constraint": None
        }
    },
    "0.1.2": {
        "model_type": "mlp",
        "input": {
            "genome": "Homo_sapiens.GRCh38.102",
            "genes": ["biotype", "protein_coding"],
        },
        "output": {
            "cl": "v2021-02-01",
            "targets": None,
        },
        "hyper_parameters": {
            "units": [256, 128],
            "activation": "selu",
            "use_bias": True,
            "l1_coef": 0.,
            "l2_coef": 0.,
            "kernel_initializer": 'lecun_normal',
            "bias_initializer": 'zeros',
            "bias_regularizer": None,
            "activity_regularizer": None,
            "kernel_constraint": None,
            "bias_constraint": None
        }
    },
    "0.1.3": {
        "model_type": "mlp",
        "input": {
            "genome": "Homo_sapiens.GRCh38.102",
            "genes": ["biotype", "protein_coding"],
        },
        "output": {
            "cl": "v2021-02-01",
            "targets": None,
        },
        "hyper_parameters": {
            "units": [512, 256, 128],
            "activation": "selu",
            "use_bias": True,
            "l1_coef": 0.,
            "l2_coef": 0.,
            "kernel_initializer": 'lecun_normal',
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
    from sfaira_extension.versions.topology_versions.human.celltype import CELLTYPEMLP_TOPOLOGIES as CELLTYPEMLP_TOPOLOGIES_EXTENSION
    CELLTYPEMLP_TOPOLOGIES = {
        **CELLTYPEMLP_TOPOLOGIES,
        **CELLTYPEMLP_TOPOLOGIES_EXTENSION
    }
except ImportError:
    pass
