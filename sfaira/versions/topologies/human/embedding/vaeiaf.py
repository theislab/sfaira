VAEIAF_TOPOLOGIES = {
    "0.1": {
        "model_type": "vaeiaf",
        "input": {
            "genome": "Homo_sapiens.GRCh38.102",
            "genes": ["biotype", "protein_coding"],
        },
        "output": {},
        "hyper_parameters": {
            "latent_dim": (256, 128, 64, 128, 256),
            "n_iaf": 2,
            "l1_coef": 0.,
            "l2_coef": 0.,
            "dropout_rate": 0.,
            "batchnorm": True,
            "activation": "tanh",
            "init": "glorot_uniform",
            "output_layer": "nb_shared_disp"
        }
    },
    "0.2": {
        "model_type": "vaeiaf",
        "input": {
            "genome": "Homo_sapiens.GRCh38.102",
            "genes": ["biotype", "protein_coding"],
        },
        "output": {},
        "hyper_parameters": {
            "latent_dim": (512, 256, 128, 256, 512),
            "n_iaf": 2,
            "l1_coef": 0.,
            "l2_coef": 0.,
            "dropout_rate": 0.,
            "batchnorm": True,
            "activation": "tanh",
            "init": "glorot_uniform",
            "output_layer": "nb_shared_disp"
        }
    }
}

# Load versions from extension if available:
try:
    from sfaira_extension.versions.topology_versions.human.embedding import VAEIAF_TOPOLOGIES as VAEIAF_TOPOLOGIES_EXTENSION
    VAEIAF_TOPOLOGIES = {
        **VAEIAF_TOPOLOGIES,
        **VAEIAF_TOPOLOGIES_EXTENSION
    }
except ImportError:
    pass
