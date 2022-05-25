import random
from typing import Callable


class _ShuffleBuffer:

    def __init__(self, generator: Callable[[], iter], buffer_size: int):
        if buffer_size < 1:
            raise ValueError('buffer_size should be larger than 0')
        self._g = generator
        self._buffer_size = buffer_size

    @staticmethod
    def buffer_replace(buffer, x):
        idx = random.randint(0, len(buffer) - 1)
        val = buffer[idx]
        buffer[idx] = x
        return val

    def iterator(self):
        buffer = []
        for x in self._g():
            if len(buffer) == self._buffer_size:
                yield self.buffer_replace(buffer, x)
            else:
                buffer.append(x)

        random.shuffle(buffer)
        while buffer:
            yield buffer.pop()


class _DatasetIteratorRepeater:

    def __init__(self, generator: Callable[[], iter], n_repeats: int):
        self._g = generator
        self._n_repeats = n_repeats

    def iterator(self):
        keep_repeating = True
        n_repetitions = 0

        while keep_repeating:
            for elem in self._g():
                yield elem

            n_repetitions += 1
            keep_repeating = (n_repetitions < self._n_repeats) or (self._n_repeats <= 0)
