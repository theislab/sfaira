from .external import SuperGenomeContainer

from . import human
from . import mouse


class Topologies:

    def __init__(
            self,
            species: str,
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
                    "vae": human.embedding.vae.VAE_TOPOLOGIES,
                    "vaeiaf": human.embedding.vaeiaf.VAEIAF_TOPOLOGIES,
                    "vaevamp": human.embedding.vaevamp.VAEVAMP_TOPOLOGIES
                }
            }
        }
        self.species = species
        self.model_class = model_class
        self.model_type = model_type
        self.topology_id = topology_id
        assert species in list(self.topologies.keys()), \
            "species %s not found in %s" % \
            (species, list(self.topologies.keys()))
        assert model_class in list(self.topologies[species].keys()), \
            "model_class %s not found in %s" % \
            (model_type, list(self.topologies[species].keys()))
        assert model_type in list(self.topologies[species][model_class].keys()), \
            "model_type %s not found in %s" % \
            (model_type, list(self.topologies[species][model_class].keys()))
        assert topology_id in list(self.topologies[species][model_class][model_type].keys()), \
            "topology_id %s not found in %s" % \
            (topology_id, list(self.topologies[species][model_class][model_type].keys()))
        self.genome_container = SuperGenomeContainer(species=species, genome=self.topology["genome"])

    @property
    def topology(self):
        return self.topologies[self.species][self.model_class][self.model_type][self.topology_id]

    @property
    def ngenes(self):
        return self.genome_container.ngenes
