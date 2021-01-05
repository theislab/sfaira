import sfaira
import sys
import tensorflow as tf

print(tf.__version__)

# Set global variables.
print("sys.argv", sys.argv)

path = str(sys.argv[1])
path_meta = str(sys.argv[2])
processes = int(float(str(sys.argv[2])))

parallelise_across_studies = True

ds = sfaira.data.dataloaders.DatasetSuperGroupSfaira(
    path=path, meta_path=path_meta, cache_path=path_meta
)
if parallelise_across_studies and processes > 1:
    ds.load_all(
        celltype_version=None,
        annotated_only=False,
        match_to_reference=None,
        remove_gene_version=True,
        load_raw=True,
        allow_caching=False,
        processes=processes,
    )
for dsg in list(ds.dataset_groups):
    if not parallelise_across_studies and processes > 1:
            # Need explicit call to load_all to uses multiprocess, otherwise this would be automatically  evoced in meta data
            # writing. Note: paralellisation is once per data set group (ie study, so that not all data is loaded into
            # memory!
            ds.load_all(
                celltype_version=None,
                annotated_only=False,
                match_to_reference=None,
                remove_gene_version=True,
                load_raw=True,
                allow_caching=False,
                processes=processes,
            )
    for k in dsg.ids:
        print(k)
        dsg.datasets[k].write_meta(dir_out=path_meta)
        # Free memory:
        del dsg.datasets[k]
