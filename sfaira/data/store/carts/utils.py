import pandas as pd


def split_batch(x):
    """
    Splits retrieval batch into consumption batches of length 1.

    Often, end-user consumption batches would be observation-wise, ie yield a first dimension of length 1.

    :param x: Data tuple of length 1 or 2: (input,) or (input, output,), where both input and output are also
         a tuple, but of batch-dimensioned tensors.
    """
    batch_dim = x[0][0].shape[0]
    for i in range(batch_dim):
        output = []
        for y in x:
            if isinstance(y, tuple):
                output.append(tuple([z.iloc[[i], :] if isinstance(z, pd.DataFrame) else z[i, :] for z in y]))
            else:
                output.append(y.iloc[[i], :] if isinstance(y, pd.DataFrame) else y[i, :])
        yield tuple(output)
