#
import numpy as np

#
import matplotlib.ticker as ticker

#
import styles


def ticks_format(value, index, interval = 2):
    """
    get the value and returns the exponent relative to the floor integer.
    """
    exp = np.log10(value)
    label = 10.**(exp%1.)

    if round(label)%interval == 0.:
        return '%i'%round(label)
    else:
        return ''


class FuncFormatter2(ticker.Formatter):
    """
    Use a user-defined function for formatting.
    The function should take in two inputs (a tick value ``x`` and a
    position ``pos``), and return a string containing the corresponding
    tick label.
    """
    def __init__(self, func, **kwargs):
        self.func = func

        self.kwargs = kwargs

    def __call__(self, x, pos=None):
        """
        Return the value of the user defined function.
        `x` and `pos` are passed through as-is.
        """
        return self.func(x, pos, **self.kwargs)
