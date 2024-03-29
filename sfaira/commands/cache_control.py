import shutil
import time

from sfaira import settings
from sfaira.consts import AdataIdsSfaira, OC


class CacheControl:

    def __init__(self):
        self.adata_ids = AdataIdsSfaira()

    @staticmethod
    def clear(countdown=True):
        if countdown:
            for i in range(5):
                print(f"WARNING: deleting entire folder {settings.cachedir_base} in {5 - i} sec ...")
                time.sleep(1)
        shutil.rmtree(path=settings.cachedir_base)

    def reload(self):
        for x in self.adata_ids.ontology_constrained:
            OC.reload_ontology(attr=x)
