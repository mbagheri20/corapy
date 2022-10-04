import numpy as np

def lorentzian(w, wj, delta=4.0):
    """Broadening by Lorentzian function

    Parameters:

    w: array
     list of frequencies
    wj: array
     frequency of centered delta function
    delta: float
     peak width, default= 4 cm^-1

    Returns:

    broadened by Lorentzian function
    """
    return 1./np.pi * (delta/2)/ ( (w - wj)**2 + (delta/2)**2 )

def gaussian(w, wj, delta=4.0):
    """Broadening by Gaussian function

    Parameters:

    w: array
     list of frequencies
    wj: array
     frequency of centered delta function
    delta: float
     peak width, default= 4 cm^-1

    Returns:

    broadened by Gaussian function
    """
    c = delta / (2 * np.sqrt(2 * np.log(2)))
    return np.exp(-np.power(w - wj, 2) / (2 * c ** 2))
