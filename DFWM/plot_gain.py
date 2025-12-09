import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import gnlse

mpl.rcParams['font.family'] = 'Times New Roman'
mpl.rcParams['font.serif'] = ['Times New Roman']
mpl.rcParams['axes.unicode_minus'] = False

def plot_gain(solution_):
    _cmap_1d = "viridis"
    _cmap_2d = "seismic"

    # plt.figure(figsize=(4, 4), facecolor='w', edgecolor='k')
    gnlse.plot_wavelength_vs_distance_logarithmic(solution_, WL_range=[300, 4000], cmap = _cmap_2d,
                                                  plot_gain=True, use_zero_norm=True,
                                                  norm_scale_inverse_factor=10)
    plt.tight_layout()
    plt.show()

    # plt.figure(figsize=(4, 4), facecolor='w', edgecolor='k')
    gnlse.plot_wavelength_vs_distance_logarithmic(solution_, WL_range=[1000, 1100], cmap = _cmap_2d,
                                                  plot_gain=True, use_zero_norm=True,
                                                  norm_scale_inverse_factor=10)
    plt.tight_layout()
    plt.show()

    # _, ax = plt.subplots()
    # gnlse.plot_wavelength_for_distance_slice_logarithmic(solution_, WL_range=[3400, 3700], ax = ax, norm = None,
    #                                                      cmap="Blues", plot_gain=True)
    # ax.set_ylim([-100, 5])
    # ax.grid(True)

    _, ax = plt.subplots()
    gnlse.plot_wavelength_for_distance_slice_logarithmic(solution_, WL_range=[3400, 3700], ax = ax,
                                                         z_slice = solution_.Z[0::10], norm = None, cmap=_cmap_1d,
                                                         plot_gain=True)
    # ax.set_ylim([-200, 50])
    ax.grid(True)

    _, ax = plt.subplots()
    gnlse.plot_wavelength_for_distance_slice_logarithmic(solution_, WL_range=[620, 640], ax = ax,
                                                         z_slice = solution_.Z[0::10], norm = None, cmap=_cmap_1d,
                                                         plot_gain=True)
    # ax.set_ylim([-200, 50])
    ax.grid(True)

    _, ax = plt.subplots()
    gnlse.plot_wavelength_for_distance_slice_logarithmic(solution_, WL_range=[1050, 1080], ax = ax,
                                                         z_slice = solution_.Z[0::10], norm = None, cmap=_cmap_1d,
                                                         plot_gain=True)
    ax.set_ylim([-800, 1400])
    ax.grid(True)
    plt.show()

    _, ax = plt.subplots()
    gnlse.visualization.plot_wavelength_slice_vs_distance_logarithmic(solution_,
                                                                      wavelengths=[625, 1064, 3550],
                                                                      ax = ax,
                                                                      norm = None, cmap=_cmap_1d, plot_gain=True)
    ax.set_xlim([0, max(solution_.Z)])
    ax.grid(True)
    plt.show()


if __name__ == "__main__":
    path = os.path.join("solutions",
        "far_detuned_fwm" + \
        "_resolution_14" + \
        "_time_window_10" + \
        "_fiber_length_0.15" +\
        "_samples_1000"+\
        ".json")

    solution = gnlse.Solution()
    solution.from_file(path)
    plot_gain(solution)