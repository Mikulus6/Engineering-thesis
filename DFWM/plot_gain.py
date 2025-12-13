import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import gnlse
import math

mpl.rcParams['font.family'] = 'Times New Roman'
mpl.rcParams['font.serif'] = ['Times New Roman']
mpl.rcParams['axes.unicode_minus'] = False

def plot_gain(solution_):
    _cmap_1d = "viridis"
    _cmap_2d = "seismic"

    plt.figure(figsize=(8, 5), facecolor='w', edgecolor='k')
    gnlse.plot_wavelength_vs_distance_logarithmic(solution_, WL_range=[300, 4000], cmap = _cmap_2d,
                                                  plot_gain=True, use_zero_norm=True,
                                                  norm_scale_inverse_factor=10)
    plt.tight_layout()
    plt.show()

    plt.figure(figsize=(8, 5), facecolor='w', edgecolor='k')
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

    fig = plt.figure(figsize=(8, 5), constrained_layout=True)
    ax = fig.add_subplot(111)
    gnlse.plot_wavelength_for_distance_slice_logarithmic(solution_, WL_range=[3400, 3700], ax = ax,
                                                         z_slice = [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4],
                                                         norm = None, cmap=_cmap_1d,
                                                         plot_gain=True)
    ax.set_ylim([-750, 2500])
    ax.legend(loc='upper right')
    ax.grid(True)

    fig = plt.figure(figsize=(8, 5), constrained_layout=True)
    ax = fig.add_subplot(111)
    gnlse.plot_wavelength_for_distance_slice_logarithmic(solution_, WL_range=[620, 640], ax = ax,
                                                         z_slice = [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4],
                                                         norm = None, cmap=_cmap_1d,
                                                         plot_gain=True)
    ax.set_ylim([-1000, 3500])
    ax.grid(True)

    fig = plt.figure(figsize=(8, 5), constrained_layout=True)
    ax = fig.add_subplot(111)
    gnlse.plot_wavelength_for_distance_slice_logarithmic(solution_, WL_range=[1050, 1080], ax = ax,
                                                        z_slice = [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4],
                                                         norm = None, cmap=_cmap_1d,
                                                         plot_gain=True)
    ax.set_ylim([-500, 2000])
    ax.grid(True)

    plt.legend(loc="upper right")
    plt.show()

    # ------------------------------------

    fig = plt.figure(figsize=(8, 5), constrained_layout=True)
    ax = fig.add_subplot(111)
    gnlse.plot_wavelength_for_distance_slice_logarithmic(solution_, WL_range=[3400, 3700], ax = ax,
                                                         z_slice = [0, 0.03, 0.06, 0.09, 0.12, 0.15],
                                                         norm = None, cmap=_cmap_1d,
                                                         plot_gain=True)
    ax.grid(True)

    fig = plt.figure(figsize=(8, 5), constrained_layout=True)
    ax = fig.add_subplot(111)
    gnlse.plot_wavelength_for_distance_slice_logarithmic(solution_, WL_range=[620, 640], ax = ax,
                                                         z_slice = [0, 0.03, 0.06, 0.09, 0.12, 0.15],
                                                         norm = None, cmap=_cmap_1d,
                                                         plot_gain=True)
    ax.grid(True)

    fig = plt.figure(figsize=(8, 5), constrained_layout=True)
    ax = fig.add_subplot(111)
    gnlse.plot_wavelength_for_distance_slice_logarithmic(solution_, WL_range=[1050, 1080], ax = ax,
                                                        z_slice = [0, 0.03, 0.06, 0.09, 0.12, 0.15],
                                                         norm = None, cmap=_cmap_1d,
                                                         plot_gain=True)
    ax.set_ylim([-1250, 1750])
    ax.grid(True)

    plt.legend(loc="upper right")
    plt.show()

    # ------------------------------------


    fig = plt.figure(figsize=(8, 5), constrained_layout=True)
    ax = fig.add_subplot(111)

    # wyidealizowana funkcja z dopasowania do teorii (porównanie)
    b = 645.45
    c = -120
    x_0 = 0.2815
    a = 0.024
    eye_guide = lambda x: (b-c)/(1+math.exp(4*(x-x_0)/a)) + c
    ax.plot(solution_.Z, [eye_guide(x) for x in solution_.Z], 'k--', lw=1, label="krzywa trendu")

    gnlse.visualization.plot_wavelength_slice_vs_distance_logarithmic(solution_,
                                                                      wavelengths=[626],
                                                                      ax = ax,
                                                                      norm = None, cmap="summer", plot_gain=True)

    ax.set_xlim([0.01, max(solution_.Z)])
    ax.set_ylim([-400, 750])
    ax.grid(True)
    plt.show()


if __name__ == "__main__":
    path = os.path.join("solutions",
        "far_detuned_fwm" + \
        "_resolution_14" + \
        "_time_window_10" + \
        "_fiber_length_0.4" +\
        "_samples_1000"+\
        ".json")

    solution = gnlse.Solution()
    solution.from_file(path)
    plot_gain(solution)