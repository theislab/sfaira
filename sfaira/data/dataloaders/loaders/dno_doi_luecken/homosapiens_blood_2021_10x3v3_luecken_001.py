import anndata
import gzip
import os
import shutil
from tempfile import TemporaryDirectory


def load(data_dir, sample_fn, **kwargs):
    fn = os.path.join(data_dir, sample_fn)
    with TemporaryDirectory() as tmpdir:
        tmppth = tmpdir + "/decompressed.h5ad"
        with gzip.open(fn, "rb") as input_f, open(tmppth, "wb") as output_f:
            shutil.copyfileobj(input_f, output_f)
        adata = anndata.read_h5ad(tmppth)
    adata.var["feature_types"] = [
        {"ATAC": "peak", "GEX": "rna", "ADT": "protein"}[x]
        for x in adata.var["feature_types"].values
    ]
    # There is NAN entries in gene_id for some antibodies:
    adata.var["gene_id"] = [str(x) for x in adata.var["gene_id"].values]
    return adata
