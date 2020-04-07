import numpy as np
from scipy import stats

def normalize_1D(array, center, scale):
    """
    Function to normalize a 1 dimensional array using the following formula:

    (array - center)/ scale

    Args:
        array - 1 dimensional numpy array
        center - a single value for the center (commonly mean(array))
        scale - a single value for the scale (commonly std(array))
    Returns:
        outarray - numpy array normalized
    """

    return (array - center)/scale

def robustz_1D(array, ignore_nan = True):
    """
    Function to Robust Z normalize an array. I.e.

    RZ = array - median(array)/(1.4826* MAD)

    Args:
        array - 1 dimensional numpy array
        ignore_nan - boolean, ignore individual nans and propogate those nans to final array
    Returns:
        outarray - numpy array robustZ normalized
    """
    if ignore_nan:
        MAD = stats.median_absolute_deviation(array, nan_policy='omit')
        median = np.nanmedian(array)
    else:
        MAD = stats.median_absolute_deviation(array)
        median = np.median(array)
    return normalize_1D(array, median, MAD)

