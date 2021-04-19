from . import human
from . import mouse
from .class_interface import TopologyContainer

TOPOLOGIES = {
    "mouse": {
        "celltype": {
            "marker": mouse.celltype.celltypemarker.CELLTYPEMARKER_TOPOLOGIES,
            "mlp": mouse.celltype.celltypemlp.CELLTYPEMLP_TOPOLOGIES
        },
        "embedding": {
            "ae": mouse.embedding.ae.AE_TOPOLOGIES,
            "linear": mouse.embedding.linear.LINEAR_TOPOLOGIES,
            "nmf": mouse.embedding.nmf.NMF_TOPOLOGIES,
            "vae": mouse.embedding.vae.VAE_TOPOLOGIES,
            "vaeiaf": mouse.embedding.vaeiaf.VAEIAF_TOPOLOGIES,
            "vaevamp": mouse.embedding.vaevamp.VAEVAMP_TOPOLOGIES
        }
    },
    "human": {
        "celltype": {
            "marker": human.celltype.celltypemarker.CELLTYPEMARKER_TOPOLOGIES,
            "mlp": human.celltype.celltypemlp.CELLTYPEMLP_TOPOLOGIES
        },
        "embedding": {
            "ae": human.embedding.ae.AE_TOPOLOGIES,
            "linear": human.embedding.linear.LINEAR_TOPOLOGIES,
            "nmf": human.embedding.nmf.NMF_TOPOLOGIES,
            "vae": human.embedding.vae.VAE_TOPOLOGIES,
            "vaeiaf": human.embedding.vaeiaf.VAEIAF_TOPOLOGIES,
            "vaevamp": human.embedding.vaevamp.VAEVAMP_TOPOLOGIES
        }
    }
}
