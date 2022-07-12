from typing import Union, Tuple, Dict

import pandas as pd


def split_batch(x: Union[Tuple, Dict]):
    """
    Splits retrieval batch into consumption batches of length 1.

    Often, end-user consumption batches would be observation-wise, ie yield a first dimension of length 1.

    :param x: Tuple or Dict
        One of the following:
            * Data tuple of length 1 or 2: (input,) or (input, output,), where both input and output are also a tuple,
            but of batch-dimensioned tensors.
            * Dict
    """
    if isinstance(x, tuple):
        batch_dim = x[0][0].shape[0]
        for i in range(batch_dim):
            output = []
            for y in x:
                if isinstance(y, tuple):
                    output.append(tuple([z.iloc[[i], :] if isinstance(z, pd.DataFrame) else z[i, :] for z in y]))
                else:
                    output.append(y.iloc[[i], :] if isinstance(y, pd.DataFrame) else y[i, :])
            yield tuple(output)
    elif isinstance(x, dict):
        keys = list(x.keys())
        batch_dim = x[keys[0]].shape[0]
        for i in range(batch_dim):
            yield {key: x[key].iloc[[i], :] if isinstance(x[key], pd.DataFrame) else x[key][i, :] for key in keys}
    else:
        raise ValueError('Input to split_batch(x) has to be either a Tuple or a Dict')
