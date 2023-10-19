from typing import Union

from sfaira.consts import OC
from sfaira.versions.metadata import OntologyHierarchical, Ontology


def get_ontology(k, organism: Union[None, str] = None) -> Union[OntologyHierarchical, None]:
    # Use global instance of ontology container:
    ocs = OC

    x = getattr(ocs, k) if hasattr(ocs, k) else None
    if x is not None and isinstance(x, dict):
        assert isinstance(organism, str), organism
        # Check if organism-specific option is available, otherwise choose generic option:
        if organism in x.keys():
            k = organism
        else:
            k = ocs.key_other
            assert k in x.keys(), x.keys()  # Sanity check on dictionary keys.
        x = x[k]
        assert x is None or isinstance(x, Ontology), x  # Sanity check on dictionary element.
    return x
