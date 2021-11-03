from . import homosapiens
from . import musmusculus
from .class_interface import TopologyContainer

TOPOLOGIES = {
    "Mus musculus": {
        "celltype": {
            "marker": musmusculus.celltype.celltypemarker.CELLTYPEMARKER_TOPOLOGIES,
            "mlp": musmusculus.celltype.celltypemlp.CELLTYPEMLP_TOPOLOGIES
        },
        "embedding": {
            "ae": musmusculus.embedding.ae.AE_TOPOLOGIES,
            "linear": musmusculus.embedding.linear.LINEAR_TOPOLOGIES,
            "nmf": musmusculus.embedding.nmf.NMF_TOPOLOGIES,
            "vae": musmusculus.embedding.vae.VAE_TOPOLOGIES,
            "vaeiaf": musmusculus.embedding.vaeiaf.VAEIAF_TOPOLOGIES,
            "vaevamp": musmusculus.embedding.vaevamp.VAEVAMP_TOPOLOGIES
        }
    },
    "Homo sapiens": {
        "celltype": {
            "marker": homosapiens.celltype.celltypemarker.CELLTYPEMARKER_TOPOLOGIES,
            "mlp": homosapiens.celltype.celltypemlp.CELLTYPEMLP_TOPOLOGIES
        },
        "embedding": {
            "ae": homosapiens.embedding.ae.AE_TOPOLOGIES,
            "linear": homosapiens.embedding.linear.LINEAR_TOPOLOGIES,
            "nmf": homosapiens.embedding.nmf.NMF_TOPOLOGIES,
            "vae": homosapiens.embedding.vae.VAE_TOPOLOGIES,
            "vaeiaf": homosapiens.embedding.vaeiaf.VAEIAF_TOPOLOGIES,
            "vaevamp": homosapiens.embedding.vaevamp.VAEVAMP_TOPOLOGIES
        }
    }
}
