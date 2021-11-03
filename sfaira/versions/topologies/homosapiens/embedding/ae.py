AE_TOPOLOGIES = {
    "0.1": {
        "model_type": "ae",
        "input": {
            "genome": "Homo_sapiens.GRCh38.102",
            "genes": ["biotype", "protein_coding"],
        },
        "output": {},
        "hyper_parameters": {
            "latent_dim": (512, 64, 512),
            "l1_coef": 0.,
            "l2_coef": 0.,
            "dropout_rate": 0.,
            "input_dropout": 0.,
            "batchnorm": True,
            "activation": "selu",
            "init": "lecun_normal",
            "output_layer": "nb_shared_disp"
        }
    },

    "0.2": {
        "model_type": "ae",
        "input": {
            "genome": "Homo_sapiens.GRCh38.102",
            "genes": ["biotype", "protein_coding"],
        },
        "output": {},
        "hyper_parameters": {
            "latent_dim": (256, 128, 64, 128, 256),
            "l1_coef": 0.,
            "l2_coef": 0.,
            "dropout_rate": 0.,
            "input_dropout": 0.,
            "batchnorm": True,
            "activation": "selu",
            "init": "lecun_normal",
            "output_layer": "nb_shared_disp"
        }
    },

    "0.3": {
        "model_type": "ae",
        "input": {
            "genome": "Homo_sapiens.GRCh38.102",
            "genes": ["biotype", "protein_coding"],
        },
        "output": {},
        "hyper_parameters": {
            "latent_dim": (512, 256, 128, 256, 512),
            "l1_coef": 0.,
            "l2_coef": 0.,
            "dropout_rate": 0.,
            "input_dropout": 0.,
            "batchnorm": True,
            "activation": "selu",
            "init": "lecun_normal",
            "output_layer": "nb_shared_disp"
        }
    },

    "0.4": {
        "model_type": "ae",
        "input": {
            "genome": "Homo_sapiens.GRCh38.102",
            "genes": ["biotype", "protein_coding"],
        },
        "output": {},
        "hyper_parameters": {
            "latent_dim": (512, 256, 128, 64, 128, 256, 512),
            "l2_coef": 0.,
            "l1_coef": 0.,
            "dropout_rate": 0.,
            "input_dropout": 0.,
            "batchnorm": True,
            "activation": "selu",
            "init": "lecun_normal",
            "output_layer": "nb_const_disp"
        }
    }
}

# Load versions from extension if available:
try:
    from sfaira_extension.versions.topology_versions.human.embedding import AE_TOPOLOGIES as AE_TOPOLOGIES_EXTENSION
    AE_TOPOLOGIES = {
        **AE_TOPOLOGIES,
        **AE_TOPOLOGIES_EXTENSION
    }
except ImportError:
    pass
