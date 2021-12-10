import re


def doi_lint(x: str) -> bool:
    """
    Check if string matches DOI format.

    Allowed formats are:
        - an actual DOI startin with "10*"
        - a string starting with "no_doi*" indicating that a dataset without corresponding DOI is used
    """
    return re.match(r'\b10\.\d+/[\w.]+\b', x) or x.startswith("no_doi")
