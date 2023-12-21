import anndata
import os
import gzip
from tempfile import TemporaryDirectory
import shutil

# the data_dir argument will be automatically set by sfaira to the folder where your datafiles lie
# the sample_fn argument will be automatically set by sfaira to each of the sample_fns provided in the yaml top section
def load(data_dir, sample_fn, fn=None, **kwargs) -> anndata.AnnData:
    fn = os.path.join(data_dir, f'GSE198623_{sample_fn}.h5ad.gz')
    uncompressed_file_type = "h5ad"
    with TemporaryDirectory() as tmpdir:
        tmppth = tmpdir + f"/decompressed.{uncompressed_file_type}"
        with gzip.open(fn, "rb") as input_f, open(tmppth, "wb") as output_f:
            shutil.copyfileobj(input_f, output_f)
        adata = anndata.read(tmppth)
    return adata