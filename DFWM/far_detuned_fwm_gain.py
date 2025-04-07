import numpy as np
import matplotlib.pyplot as plt
import gnlse
from gnlse.common import c
from scipy.interpolate import RectBivariateSpline


def get_point(solution_: gnlse.Solution, distance, wavelength):
    distance_index = np.abs(solution_.Z - distance).argmin()
    wavelength_index = np.abs(solution_.W - wavelength).argmin()
    return solution_.AW[distance_index, wavelength_index]

if __name__ == "__main__":
    path = "far_detuned_fwm" + \
        "_resolution_14" + \
        "_time_window_20" + \
        "_fiber_length_0.1" + \
        ".json"

    solution = gnlse.Solution()
    solution.from_file(path)

    WL_range = [300, 2000]

    norm = np.max(np.abs(solution.AW)**2)

    IW = np.fliplr(np.abs(solution.AW)**2 / norm)
    WL = 2 * np.pi * c / solution.W  # wavelength grid
    WL_asc = np.flip(WL, )  # ascending order for interpolation
    iis = np.logical_and(WL_asc > WL_range[0],
                         WL_asc < WL_range[1])  # indices of interest

    WL_asc = WL_asc[iis]
    IW = IW[:, iis]

    interpolator = RectBivariateSpline(solution.Z, WL_asc, IW)
    newWL = np.linspace(np.min(WL_asc), np.max(WL_asc), IW.shape[1])
    toshow = interpolator(solution.Z, newWL)

    for i in np.linspace(0, np.max(solution.Z), 10):
        row_index = np.abs(solution.Z - i).argmin()
        row_data = IW[row_index, :]  # 1D array

        # Plot the row as a 1D line plot
        plt.plot(row_data)
    plt.grid(True)
    plt.show()
