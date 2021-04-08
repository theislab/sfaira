LINEAR_TOPOLOGIES = {
    "0.1": {
        "genome": "Mus_musculus.GRCm38.102",
        "genes": ["protein_coding"],
        "hyper_parameters": {
            "latent_dim": 64,
            "l1_coef": 0.,
            "l2_coef": 0.,
            "positive_components": False,
            "output_layer": "nb_shared_disp"
        }
    },

    "0.2": {
        "genome": "Mus_musculus.GRCm38.102",
        "genes": ["protein_coding"],
        "hyper_parameters": {
            "latent_dim": 128,
            "l1_coef": 0.,
            "l2_coef": 0.,
            "positive_components": False,
            "output_layer": "nb_shared_disp"
        }
    },

    "0.3": {
        "genome": "Mus_musculus.GRCm38.102",
        "genes": ["protein_coding"],
        "hyper_parameters": {
            "latent_dim": 128,
            "l1_coef": 0.,
            "l2_coef": 0.,
            "positive_components": False,
            "output_layer": "nb_const_disp"
        }
    }
}

# Load versions from extension if available:
try:
    from sfaira_extension.versions.topology_versions.mouse.embedding import LINEAR_TOPOLOGIES as LINEAR_TOPOLOGIES_EXTENSION
    LINEAR_TOPOLOGIES = {
        **LINEAR_TOPOLOGIES,
        **LINEAR_TOPOLOGIES_EXTENSION
    }
except ImportError:
    pass