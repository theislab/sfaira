import shutil
import time

from sfaira import settings
from sfaira.consts import AdataIdsSfaira, OCS


class CacheControl:

    def __init__(self):
        self.adata_ids = AdataIdsSfaira()

    @staticmethod
    def clear():
        for i in range(5):
            print(f"WARNING: deleting entire folder {settings.cachedir_ontologies} in {5 - i} sec ...")
            time.sleep(1)
        shutil.rmtree(path=settings.cachedir_ontologies)

    def reload(self):
        for x in self.adata_ids.ontology_constrained:
            OCS.reload_ontology(attr=x)
