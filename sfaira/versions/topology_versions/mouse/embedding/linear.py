LINEAR_TOPOLOGIES = {
    "0.1": {
        "genome": "Mus_musculus_GRCm38_97",
        "hyper_parameters": {
                "latent_dim": 64,
                "l1_coef": 0.,
                "l2_coef": 0.,
                "positive_components": False,
                "output_layer": "nb_shared_disp"
        }
    },

    "0.2": {
        "genome": "Mus_musculus_GRCm38_97",
        "hyper_parameters": {
                "latent_dim": 128,
                "l1_coef": 0.,
                "l2_coef": 0.,
                "positive_components": False,
                "output_layer": "nb_shared_disp"
        }
    },

    "0.3": {
        "genome": "Mus_musculus_GRCm38_97",
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
    import sfaira_extension.api as sfairae
    ADD_TOPOLOGIES = sfairae.versions.topology_versions.mouse.embedding.LINEAR_TOPOLOGIES
    for k in LINEAR_TOPOLOGIES.keys():
        if k in ADD_TOPOLOGIES.keys():
            LINEAR_TOPOLOGIES.update(ADD_TOPOLOGIES)
except ImportError:
    pass
