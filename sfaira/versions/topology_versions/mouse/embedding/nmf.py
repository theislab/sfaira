NMF_TOPOLOGIES = {
    "0.1": {
        "genome": "Mus_musculus_GRCm38_97",
        "hyper_parameters": {
                "latent_dim": 64,
                "l1_coef": 0.,
                "l2_coef": 0.,
                "positive_components": True,
                "output_layer": "nb_shared_disp"
        }
    },

    "0.2": {
        "genome": "Mus_musculus_GRCm38_97",
        "hyper_parameters": {
                "latent_dim": 128,
                "l1_coef": 0.,
                "l2_coef": 0.,
                "positive_components": True,
                "output_layer": "nb_shared_disp"
        }
    },

    "0.3": {
        "genome": "Mus_musculus_GRCm38_97",
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
    import sfaira_extension
    NMF_TOPOLOGIES = {
        **NMF_TOPOLOGIES,
        **sfaira_extension.versions.topology_versions.mouse.embedding.NMF_TOPOLOGIES
    }
except ImportError:
    pass
