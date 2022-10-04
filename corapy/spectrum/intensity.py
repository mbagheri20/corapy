import numpy as np
from scipy import constants


def calculate_intensity(freqs, tensors,num_activity):
    """Calculate Raman intensity and activity

    Parameters:

    freqs: array
     list of frequncies
    tensors: array
     list of raman tensors
    num_activity: Integer
     number of activities

    Returns:

    activity: array
     Raman activity
    intensity: array
     Raman intensity
    """
    alpha2 = []
    beta2 = []
    for i in range(num_activity):
        t = tensors[i]
        alpha2.append(((t[0][0] + t[1][1] + t[2][2]) / 3.0) ** 2)
        beta2.append(((t[0][0] - t[1][1]) ** 2 + (t[0][0] - t[2][2]) ** 2 + (t[1][1] - t[2][2]) ** 2 + 6.0 * (
                t[0][1] ** 2 + t[0][2] ** 2 + t[1][2] ** 2)) / 2.0)

    hbar = constants.hbar
    kT = 4.11 * 1e-21  # J
    wn = [element * 0.02998 * constants.tera * 2 * np.pi for element in freqs]  # Hz

    n = []
    for i in range(len(wn)):
        n.append(1 / (np.exp(hbar * wn[i] / kT) - 1))

    activity = []
    intensity = []

    for i in range(len(alpha2)):
        activity.append((45 * alpha2[i]) + (7 * beta2[i]))
        intensity.append(((45 * alpha2[i]) + (7 * beta2[i])) * (((n[i] + 1)) / wn[i]))

    return activity, intensity