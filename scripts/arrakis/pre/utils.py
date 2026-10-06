import numpy as np

def smooth_data(array, span):
    """
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Smooth data with moving average
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    @param array (array): input data
    @param span (int): smoothing span (should be odd)
    @return output (array): smoothed data
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    """
    pad_width = span // 2
    array_padded = np.pad(array, pad_width, mode='edge')
    output = np.convolve(array_padded, np.ones(span)/span, mode='valid')
    output = output[:len(array)]  # in case of even span
    return output
