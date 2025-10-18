import os
import matplotlib.pyplot as plt
import gnlse


def plot_solution(solution_):
    _cmap_1d = "gnuplot"
    _cmap_2d = "CMRmap"
    plt.figure(figsize=(8, 4), facecolor='w', edgecolor='k')
    plt.subplot(1, 2, 1)
    gnlse.plot_wavelength_vs_distance_logarithmic(solution_, WL_range=[300, 4000], cmap = _cmap_2d)

    plt.subplot(1, 2, 2)
    gnlse.plot_delay_vs_distance_logarithmic(solution_, time_range=[min(solution_.t),
                                                                    max(solution_.t)], cmap = _cmap_2d)

    plt.tight_layout()
    plt.show()

    plt.figure(figsize=(8, 4), facecolor='w', edgecolor='k')
    plt.subplot(1, 2, 1)
    gnlse.plot_wavelength_vs_distance_logarithmic(solution_, WL_range=[1000, 1100], cmap = _cmap_2d)

    plt.subplot(1, 2, 2)
    gnlse.plot_delay_vs_distance_logarithmic(solution_, time_range=[min(solution_.t),
                                                                    max(solution_.t)], cmap = _cmap_2d)

    plt.tight_layout()
    plt.show()

    _, ax = plt.subplots()
    gnlse.plot_wavelength_for_distance_slice_logarithmic(solution_, WL_range=[3400, 3700], ax = ax, norm = 51e3,
                                                         cmap=_cmap_1d)
    ax.set_ylim([-100, 5])
    ax.grid(True)

    _, ax = plt.subplots()
    gnlse.plot_wavelength_for_distance_slice_logarithmic(solution_, WL_range=[3400, 3700], ax = ax,
                                                         z_slice = solution_.Z[0::10], norm = 51e3, cmap=_cmap_1d)
    ax.set_ylim([-200, 50])
    ax.grid(True)

    _, ax = plt.subplots()
    gnlse.plot_wavelength_for_distance_slice_logarithmic(solution_, WL_range=[620, 640], ax = ax,
                                                         z_slice = solution_.Z[0::10], norm = 51e3, cmap=_cmap_1d)
    ax.set_ylim([-200, 50])
    ax.grid(True)

    _, ax = plt.subplots()
    gnlse.plot_wavelength_for_distance_slice_logarithmic(solution_, WL_range=[1050, 1080], ax = ax,
                                                         z_slice = solution_.Z[0::10], norm = 51e3, cmap=_cmap_1d)
    ax.set_ylim([-200, 50])
    ax.grid(True)
    plt.show()


if __name__ == "__main__":
    path = os.path.join("solutions",
        "far_detuned_fwm" + \
        "_resolution_14" + \
        "_time_window_10" + \
        "_fiber_length_0.2" + \
        ".json")

    solution = gnlse.Solution()
    solution.from_file(path)
    plot_solution(solution)