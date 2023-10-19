import os
import shutil
import tarfile
import tempfile

import scanpy as sc


def load(data_dir, sample_fn, **kwargs):
    fn = os.path.join(data_dir, 'GSE122960_RAW.tar')

    with tarfile.open(fn) as tar:
        with tempfile.TemporaryDirectory() as tmpdir:
            f = tar.extractfile(f'{sample_fn}_filtered_gene_bc_matrices_h5.h5')
            tmppth = tmpdir + "/decompressed.h5"
            with open(tmppth, "wb") as output_f:
                shutil.copyfileobj(f, output_f)
            adata = sc.read_10x_h5(tmppth)

    return adata
