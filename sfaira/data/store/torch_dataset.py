import torch


class SfairaIterableDataset(torch.utils.data.IterableDataset):

    def __init__(self, iterator_fun, **kwargs):
        super(SfairaIterableDataset, self).__init__(**kwargs)
        self.iterator_fun = iterator_fun

    def __iter__(self):
        return self.iterator_fun()
