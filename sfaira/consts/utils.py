import os


def clean_doi(doi: str):
    return f'd{doi.translate({ord(c): "_" for c in r"!@#$%^&*()[]/{};:,.<>?|`~-=_+"})}'


def clean_id_str(s):
    if s is not None:
        s = s.replace(',', '').replace(' ', '').replace('-', '').replace('_', '').replace("'", '').lower()
    return s
