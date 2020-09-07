AE_TOPOLOGIES = {
    "0.1": {
        "genome": "Homo_sapiens_GRCh38_97",
        "hyper_parameters": {
                "latent_dim": (512, 64, 512),
                "l2_coef": 0.,
                "l1_coef": 0.,
                "dropout_rate": 0.,
                "input_dropout": 0.,
                "batchnorm": False,
                "activation": "tanh",
                "init": "glorot_uniform",
                "output_layer": "nb_shared_disp"
        }
    },

    "0.2": {
        "genome": "Homo_sapiens_GRCh38_97",
        "hyper_parameters": {
                "latent_dim": (256, 128, 64, 128, 256),
                "l2_coef": 0.,
                "l1_coef": 0.,
                "dropout_rate": 0.,
                "input_dropout": 0.,
                "batchnorm": False,
                "activation": "tanh",
                "init": "glorot_uniform",
                "output_layer": "nb_shared_disp"
        }
    },

    "0.3": {
        "genome": "Homo_sapiens_GRCh38_97",
        "hyper_parameters": {
                "latent_dim": (512, 256, 128, 256, 512),
                "l2_coef": 0.,
                "l1_coef": 0.,
                "dropout_rate": 0.,
                "input_dropout": 0.,
                "batchnorm": True,
                "activation": "relu",
                "init": "glorot_uniform",
                "output_layer": "nb_shared_disp"
        }
    },

    "0.4": {
        "genome": "Homo_sapiens_GRCh38_97",
        "hyper_parameters": {
                "latent_dim": (512, 256, 128, 64, 128, 256, 512),
                "l2_coef": 0.,
                "l1_coef": 0.,
                "dropout_rate": 0.,
                "input_dropout": 0.,
                "batchnorm": True,
                "activation": "relu",
                "init": "glorot_uniform",
                "output_layer": "nb_const_disp"
        }
    }
}

# Load versions from extension if available:
try:
    import sfaira_extension as sfairae
    ADD_TOPOLOGIES = sfairae.versions.topology_versions.human.embedding.AE_TOPOLOGIES
    for k in AE_TOPOLOGIES.keys():
        if k in ADD_TOPOLOGIES.keys():
            AE_TOPOLOGIES.update(ADD_TOPOLOGIES)
except ImportError:
    pass
