VAEVAMP_TOPOLOGIES = {
    "0.1": {
        "genome": "Homo_sapiens_GRCh38_97",
        "hyper_parameters": {
            "latent_dim": (32, 32),
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
        "genome": "Homo_sapiens_GRCh38_97",
        "hyper_parameters": {
            "latent_dim": (64, 64),
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
    import sfaira_extension.api as sfairae
    ADD_TOPOLOGIES = sfairae.versions.topology_versions.human.embedding.VAEVAMP_TOPOLOGIES
    for k in VAEVAMP_TOPOLOGIES.keys():
        if k in ADD_TOPOLOGIES.keys():
            VAEVAMP_TOPOLOGIES.update(ADD_TOPOLOGIES)
except ImportError:
    pass
