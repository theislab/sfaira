import anndata
import os
import scipy.sparse


# the data_dir argument will be automatically set by sfaira to the folder where your datafiles lie
def load(data_dir, **kwargs):

    fn = os.path.join(data_dir, "GSE187877_expression_count_matrix.txt.gz")
    adata = anndata.read_csv(fn, delimiter=" ").T

    # Using Pubchem names to define treatment compounds
    treat_dict = {
        "C": "",
        "D": "Androstanolone",
        "E": "Estradiol",
    }
    state_dict = {
        "C": "",
        "D": "30 nM",
        "E": "100 nM",
    }
    sample_dict = {
        "C1": "control_1",
        "C2": "control_2",
        "D1": "dht_1",
        "D2": "dht_2",
        "E1": "oestradiol_1",
        "E2": "oestradiol_2",
    }

    treatment = []
    sample = []
    state = []
    for i in adata.obs.index:
        ii = i.split("_")
        treatment.append(treat_dict[ii[0]])
        state.append(state_dict[ii[0]])
        sample.append(sample_dict[ii[1]])

    adata.obs["treatment"] = treatment
    adata.obs["sample"] = sample
    adata.obs["treatment_concentrations"] = sample
    adata.X = scipy.sparse.csr_matrix(adata.X)

    return adata
