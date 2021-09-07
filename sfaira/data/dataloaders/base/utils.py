from typing import List, Union
from sfaira.versions.metadata import Ontology

UNS_STRING_META_IN_OBS = "__obs__"


def is_child(
        query,
        ontology: Union[Ontology, bool, int, float, str, List[bool], List[int], List[float], List[str]],
        ontology_parent=None,
) -> True:
    """
    Check whether value is from set of allowed values using ontology.

    :param query: Value to attempt to set, only yield a single value and not a list.
    :param ontology: Constraint for values.
        Either ontology instance used to constrain entries, or list of allowed values.
    :param ontology_parent: If ontology is a DAG, not only check if node is a DAG node but also whether it is a child
        of this parent node.
    :return: Whether attempted term is sub-term of allowed term in ontology
    """
    if ontology_parent is None and ontology is None:
        return True
    else:
        if isinstance(ontology, Ontology):
            if ontology_parent is None:
                return ontology.is_node(query)
            else:
                return ontology.is_a(query=query, reference=ontology_parent)
        elif ontology is None:
            return query == ontology_parent
        else:
            raise ValueError(f"did not recognize ontology type {type(ontology)}")


def clean_string(s):
    if s is not None:
        s = s.replace(',', '').replace(' ', '').replace('-', '').replace('_', '').replace("'", '').lower()
    return s


def get_directory_formatted_doi(x: str) -> str:
    return "d" + "_".join("_".join("_".join(x.split("/")).split(".")).split("-"))
