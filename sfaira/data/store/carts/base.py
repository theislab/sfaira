from functools import partial

import pandas as pd

from sfaira.data.store.batch_schedule import BatchDesignBase


class CartBase:
    """
    A cart manages data emission for a specified subset of a store.

    Data can be queried through:

        1) an iterator function
        2) the full data objects (.x and .obs)

    1) Iterator function:
    The iterator can often directly be used without a class like this one around it.
    However, that often implies that the namespace with the pointer to the data set is destroyed after iterator
    function declaration, which means that the pointer needs to be redefined for every full pass over the iterator.
    The class around this property maintains the namespace that includes the pointer and its instance can be used to
    avoid redefining the pointer every time the generator runs out and is re-called.
    For this run time advantage to not be compromised by memory load and class initialisation run time cost induced by
    actually copying data objects, it is important that the data object stored in this class is indeed a pointer.
    This is the case for:

        - lazily loaded dask arrays
        - anndata.Anndata view

    which have their own classes below.
    """

    schedule: BatchDesignBase
    var: pd.DataFrame  # Feature meta data (features x properties).
    map_fn: callable

    def adaptor(
            self,
            generator_type: str,
            dataset_kwargs: dict = None,
            shuffle_buffer: int = 0,
            repeat: int = 1,
            **kwargs
    ):
        """
        The adaptor turns a python base generator into a different iteratable object, defined by generator_type.

        :param generator_type: Type of output iteratable.
            - python base generator (no change to `.generator`)
            - tensorflow dataset: This dataset is defined on a python iterator.
                Important:
                    This only returns the tf.data.Dataset.from_generator(). You need to define the input pipeline
                    (e.g. .batch(), .prefetch()) on top of this data set.
            - pytorch: We distinguish torch.data.Dataset and torch.data.DataLoader ontop of either.
                The Dataset vs DataLoader distinction is made by the "" suffix for Dataset or "-loader" suffix for +
                dataloader. The distinction between Dataset and IteratableDataset defines if the object is defined
                directly on a dask array or based on a python iterator on a dask array. Note that the python iterator
                can implement favorable remote access schemata but the torch.data.Dataset generally causes less trouble
                in out-of-the-box usage.
                    - torch.data.Dataset: "torch" prefix, ie "torch" or "torch-loader"
                    - torch.data.IteratableDataset: "torch-iter" prefix, ie "torch-iter" or "torch-iter-loader"
                Important:
                    For model training in pytorch you need the "-loader" prefix. You can specify the arguments passed
                    to torch.utils.data.DataLoader by the dataset_kwargs dictionary.
        :param dataset_kwargs: Dict
            Parameters to pass to the constructor of torch Dataset.
            Only relevant if generator_type in ['torch', 'torch-loader']
        :param shuffle_buffer: int
            If shuffle_buffer > 0 -> Use a shuffle buffer with size shuffle_buffer to shuffle output of self.iterator
            (this option is useful when using randomized_batch_access in the DaskCart)
        :param repeat: int
            Number of times to repeat the dataset until the underlying generator runs out of samples.
            If repeat <= 0 -> repeat dataset forever
        :returns: Modified iteratable (see generator_type).
        """
        if not dataset_kwargs:
            dataset_kwargs = {}

        if generator_type == "python":
            g = self.iterator(repeat=repeat, shuffle_buffer=shuffle_buffer)
        elif generator_type == "tensorflow":
            import tensorflow as tf

            g = tf.data.Dataset.from_generator(
                generator=partial(self.iterator, repeat=repeat, shuffle_buffer=shuffle_buffer), **kwargs
            )
        elif generator_type in ["torch", "torch-loader"]:
            from torch.utils.data import DataLoader
            # Only import this module if torch is used to avoid strict torch dependency:
            from sfaira.data.store.torch_dataset import SfairaDataset

            g = SfairaDataset(map_fn=self.map_fn, obs=self.obs, x=self.x, **dataset_kwargs)
            if generator_type == "torch-loader":
                g = DataLoader(g, **kwargs)
        elif generator_type in ["torch-iter", "torch-iter-loader"]:
            from torch.utils.data import DataLoader
            # Only import this module if torch is used to avoid strict torch dependency:
            from sfaira.data.store.torch_dataset import SfairaIterableDataset

            g = SfairaIterableDataset(iterator_fun=partial(self.iterator, repeat=repeat, shuffle_buffer=shuffle_buffer))
            if generator_type == "torch-iter-loader":
                g = DataLoader(g, **kwargs)
        else:
            raise ValueError(f"{generator_type} not recognized")
        return g

    @property
    def adata(self):
        """
        Assembles a slice of this cart based on .obs_idx as an anndata instance.
        """
        raise NotImplementedError()

    def iterator(self, repeat: int = 1, shuffle_buffer: int = 0):
        """
        Iterator over data matrix and meta data table, yields batches of data points.
        """
        raise NotImplementedError()

    @property
    def obs_idx(self):
        """
        Integer observation indices to select from cart. These will be emitted if data is queried.
        """
        raise NotImplementedError()

    @property
    def n_batches(self):
        return self.schedule.n_batches

    @property
    def n_obs(self):
        """Total number of observations in cart."""
        raise NotImplementedError()

    @property
    def n_obs_selected(self):
        """Total number of selected observations in cart."""
        raise NotImplementedError()

    @property
    def n_var(self):
        """Total number of features defined for return in cart."""
        raise NotImplementedError()

    @property
    def obs(self):
        """
        Selected meta data matrix (cells x meta data) that is emitted in batches by .iterator().
        """
        raise NotImplementedError()

    def move_to_memory(self):
        """
        Load underlying array into memory into memory.
        """
        raise NotImplementedError()

    @property
    def x(self):
        """
        Selected data matrix (cells x features) that is emitted in batches by .iterator().
        """
        raise NotImplementedError()
