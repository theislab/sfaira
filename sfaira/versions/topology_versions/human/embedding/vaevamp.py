VAEVAMP_TOPOLOGIES = {
    "0.2": {
        "genome": "Homo_sapiens_GRCh38_97",
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
        "genome": "Homo_sapiens_GRCh38_97",
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
    import sfaira_extension
    ADD_TOPOLOGIES = sfaira_extension.versions.topology_versions.human.embedding.VAEVAMP_TOPOLOGIES
    for k in VAEVAMP_TOPOLOGIES.keys():
        if k in ADD_TOPOLOGIES.keys():
            VAEVAMP_TOPOLOGIES.update(ADD_TOPOLOGIES)
except ImportError:
    pass
