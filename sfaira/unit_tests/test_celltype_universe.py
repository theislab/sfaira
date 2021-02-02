import numpy as np
import pandas as pd
import unittest

from sfaira.versions.celltype_versions import OntologyObo, ORGANISM_DICT


class TestCellTypeUniverse(unittest.TestCase):
    dir_debugging = "~/Desktop/temp/"
    dir_debugging2 = "~/Desktop/temp2/"
    dir_debugging3 = "~/Desktop/temp3/"

    def test_debugging(self, reduced=False):
        import csv
        onto = OntologyObo()
        for k, v in ORGANISM_DICT.items():
            for kk, vv in v.items():
                universe = vv.celltype_universe["0"]
                tab = onto.fuzzy_match_nodes(universe, match_only=True)
                if not np.all(tab["matched"].values):
                    tab2 = onto.fuzzy_match_nodes(universe, match_only=False, include_old=True, remove=["unkown"])
                    if not reduced:
                        tab2.to_csv(
                            self.dir_debugging + k + "_" + kk + "_universe.csv",
                            index=False, quoting=csv.QUOTE_NONE, sep=";"
                        )
                    else:
                        tab2.loc[not tab["matched"].values].to_csv(
                            self.dir_debugging + k + "_" + kk + "_universe.csv",
                            index=False, quoting=csv.QUOTE_NONE
                        )

    def test_debugging2(self):
        import csv
        onto = OntologyObo()
        for k, v in ORGANISM_DICT.items():
            for kk, vv in v.items():
                names = list(vv.ontology["0"]["names"].keys())
                tab = onto.fuzzy_match_nodes(names, match_only=True)
                if not np.all(tab["matched"].values):
                    tab = onto.fuzzy_match_nodes(names, match_only=False, include_old=True, remove=["unkown"])
                    tab.to_csv(
                        self.dir_debugging2 + k + "_" + kk + "_universe.csv",
                        index=False, quoting=csv.QUOTE_NONE, sep=";"
                    )

    def test_debugging3(self):
        import csv
        onto = OntologyObo()
        tab = pd.DataFrame({"name,id": [",".join([x, y]) for x, y in zip(
            [v["name"] for k, v in onto.graph.nodes.items()],
            list(onto.graph.nodes.keys())
        )]})
        tab.to_csv(
            self.dir_debugging3 + "onto_full.csv",
            index=False, quoting=csv.QUOTE_NONE, sep=";"
        )

    def test_only(self):
        onto = OntologyObo()
        for k, v in ORGANISM_DICT.items():
            for kk, vv in v.items():
                universe = vv.celltype_universe["0"]
                tab = onto.fuzzy_match_nodes(universe, match_only=True)
                print(tab.loc[not tab["matched"].values])
                assert np.all(tab["matched"].values), f"{k} {kk}"


if __name__ == '__main__':
    unittest.main()
