import matplotlib.pyplot as plt
import numpy as np
from corapy.spectrum.functions import lorentzian, gaussian



def plot_spectrum(freqs, intensities, width, func='lorentzian', file_name=None):
    """Calculate Raman intensity and activity

    Parameters:

    freqs: array
     list of frequencies in cm^-1
    intensities: array
     Raman intensity or activity
    width: float
     width of the Gaussian or Lorentzian function
    func: string
     name of broadening function(lorentzian or gaussian), default= lorentzian
    file_name: string
     name file for saving plot with extension (e.g. plot.png)
    """
    wavenumbers = np.linspace(0, np.max(freqs) + 50, 1000)
    broadened_intensities = np.zeros_like(wavenumbers)

    for i, j in enumerate(freqs):
        if func=='lorentzian':
            broadened_intensities += intensities[i] * lorentzian(wavenumbers, j, width)
        else:
            broadened_intensities += intensities[i] * gaussian(wavenumbers, j, width)

    plt.plot(wavenumbers, np.array(broadened_intensities)/np.max(broadened_intensities), 'r')

    plt.ylabel('Raman intensity (a.u.)', fontsize=18)
    plt.xlabel('Raman shift (cm$^{-1}$)', fontsize=18)
    plt.xlim(np.min(wavenumbers), np.max(wavenumbers))
    plt.yticks(fontsize=18)
    plt.xticks(fontsize=18)
    #plt.legend(fontsize=18)
    plt.tight_layout()

    if file_name is not None:
        plt.savefig(file_name)
    else:
        plt.show()