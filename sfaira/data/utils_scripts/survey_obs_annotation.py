import numpy as np
import sfaira
import sys

# Set global variables.
print("sys.argv", sys.argv)

data_path = str(sys.argv[1])
path_meta = str(sys.argv[2])
path_cache = str(sys.argv[3])

universe = sfaira.data.dataloaders.Universe(
    data_path=data_path, meta_path=path_meta, cache_path=path_cache
)
for k, v in universe.datasets.items():
    print(k)
    v.load(
        load_raw=False,
        allow_caching=True,
    )
    for col in v.adata.obs.columns:
        val = np.sort(np.unique(v.adata.obs[col].values))
        if len(val) > 20:
            val = val[:20]
        print(f"{k}: {col}: {val}")
    v.clear()
