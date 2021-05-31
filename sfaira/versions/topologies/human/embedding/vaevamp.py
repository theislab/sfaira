VAEVAMP_TOPOLOGIES = {
    "0.2": {
        "model_type": "vaevamp",
        "input": {
            "genome": "Homo_sapiens.GRCh38.102",
            "genes": ["biotype", "protein_coding"],
        },
        "output": {},
        "hyper_parameters": {
            "latent_dim": (256, 128, (32, 32), 128, 256),
            "l1_coef": 0.,
            "l2_coef": 0.,
            "dropout_rate": 0.,
            "batchnorm": True,
            "activation": "selu",
            "init": "lecun_normal",
            "output_layer": "nb_shared_disp"
        }
    },
    "0.3": {
        "model_type": "vaevamp",
        "input": {
            "genome": "Homo_sapiens.GRCh38.102",
            "genes": ["biotype", "protein_coding"],
        },
        "output": {},
        "hyper_parameters": {
            "latent_dim": (512, 256, (64, 64), 256, 512),
            "l1_coef": 0.,
            "l2_coef": 0.,
            "dropout_rate": 0.,
            "batchnorm": True,
            "activation": "selu",
            "init": "lecun_normal",
            "output_layer": "nb_shared_disp"
        }
    }
}

# Load versions from extension if available:
try:
    from sfaira_extension.versions.topology_versions.human.embedding import VAEVAMP_TOPOLOGIES as VAEVAMP_TOPOLOGIES_EXTENSION
    VAEVAMP_TOPOLOGIES = {
        **VAEVAMP_TOPOLOGIES,
        **VAEVAMP_TOPOLOGIES_EXTENSION
    }
except ImportError:
    pass
