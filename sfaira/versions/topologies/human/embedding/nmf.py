NMF_TOPOLOGIES = {
    "0.1": {
        "model_type": "linear",
        "input": {
            "genome": "Homo_sapiens.GRCh38.102",
            "genes": ["biotype", "protein_coding"],
        },
        "output": {},
        "hyper_parameters": {
            "latent_dim": 64,
            "l1_coef": 0.,
            "l2_coef": 0.,
            "positive_components": True,
            "output_layer": "nb_shared_disp"
        }
    },

    "0.2": {
        "model_type": "linear",
        "input": {
            "genome": "Homo_sapiens.GRCh38.102",
            "genes": ["biotype", "protein_coding"],
        },
        "output": {},
        "hyper_parameters": {
            "latent_dim": 128,
            "l1_coef": 0.,
            "l2_coef": 0.,
            "positive_components": True,
            "output_layer": "nb_shared_disp"
        }
    },

    "0.3": {
        "model_type": "linear",
        "input": {
            "genome": "Homo_sapiens.GRCh38.102",
            "genes": ["biotype", "protein_coding"],
        },
        "output": {},
        "hyper_parameters": {
            "latent_dim": 128,
            "l1_coef": 0.,
            "l2_coef": 0.,
            "positive_components": True,
            "output_layer": "nb_const_disp"
        }
    }
}

# Load versions from extension if available:
try:
    from sfaira_extension.versions.topology_versions.human.embedding import NMF_TOPOLOGIES as NMF_TOPOLOGIES_EXTENSION
    NMF_TOPOLOGIES = {
        **NMF_TOPOLOGIES,
        **NMF_TOPOLOGIES_EXTENSION
    }
except ImportError:
    pass
