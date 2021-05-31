CELLTYPEMARKER_TOPOLOGIES = {
    "0.0.1": {
        "model_type": "marker",
        "input": {
            "genome": "Homo_sapiens.GRCh38.102",
            "genes": ["biotype", "protein_coding"],
        },
        "output": {
            "cl": "v2021-02-01",
            "targets": None,
        },
        "hyper_parameters": {
            "l1_coef": 0.,
            "l2_coef": 0.,
            "kernel_initializer": 'glorot_uniform',
            "bias_initializer": 'zeros',
            "bias_regularizer": None,
            "kernel_constraint": None,
            "bias_constraint": None
        }
    }
}

# Load versions from extension if available:
try:
    from sfaira_extension.versions.topology_versions.human.celltype import CELLTYPEMARKER_TOPOLOGIES as CELLTYPEMARKER_TOPOLOGIES_EXTENSION
    CELLTYPEMARKER_TOPOLOGIES = {
        **CELLTYPEMARKER_TOPOLOGIES,
        **CELLTYPEMARKER_TOPOLOGIES_EXTENSION
    }
except ImportError:
    pass
