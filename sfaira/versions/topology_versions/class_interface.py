from .external import SuperGenomeContainer

from . import human
from . import mouse


class Topologies:

    def __init__(
            self,
            organism: str,
            model_class: str,
            model_type: str,
            topology_id: str
    ):
        self.topologies = {
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
            "loaders": {
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
        self.organism = organism
        self.model_class = model_class
        self.model_type = model_type
        self.topology_id = topology_id
        assert organism in list(self.topologies.keys()), \
            "organism %s not found in %s" % \
            (organism, list(self.topologies.keys()))
        assert model_class in list(self.topologies[organism].keys()), \
            "model_class %s not found in %s" % \
            (model_type, list(self.topologies[organism].keys()))
        assert model_type in list(self.topologies[organism][model_class].keys()), \
            "model_type %s not found in %s" % \
            (model_type, list(self.topologies[organism][model_class].keys()))
        assert topology_id in list(self.topologies[organism][model_class][model_type].keys()), \
            "topology_id %s not found in %s" % \
            (topology_id, list(self.topologies[organism][model_class][model_type].keys()))
        self.genome_container = SuperGenomeContainer(organism=organism, genome=self.topology["genome"])

    @property
    def topology(self):
        return self.topologies[self.organism][self.model_class][self.model_type][self.topology_id]

    @property
    def ngenes(self):
        return self.genome_container.ngenes
