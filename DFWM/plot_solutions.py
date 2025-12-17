import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import gnlse
from DFWM.plot_gain import plot_gain
from matplotlib.ticker import FuncFormatter

mpl.rcParams['font.family'] = 'Times New Roman'
mpl.rcParams['font.serif'] = ['Times New Roman']
mpl.rcParams['axes.unicode_minus'] = False

def comma_only_formatter(axis, decimals=2):
    base = axis.get_major_formatter()
    def fmt(x, pos):
        s = base.format_data_short(x)
        return f"{float(s):.{decimals}f}".replace('.', ',') if decimals > 0 else str(int(round(float(s))))
    return FuncFormatter(fmt)

def plot_solution(solution_):
    _cmap_1d = "viridis"
    _cmap_2d = "CMRmap"
    # plt.figure(figsize=(8, 4), facecolor='w', edgecolor='k')
    # plt.subplot(1, 2, 1)

    plt.figure(figsize=(8, 5), facecolor='w', edgecolor='k')
    gnlse.plot_wavelength_vs_distance_logarithmic(solution_, WL_range=[300, 4000], cmap = _cmap_2d)

    # plt.subplot(1, 2, 2)
    # gnlse.plot_delay_vs_distance_logarithmic(solution_, time_range=[min(solution_.t),
    #                                                                 max(solution_.t)], cmap = _cmap_2d)

    ax = plt.gca()
    ax.xaxis.set_major_formatter(comma_only_formatter(ax.xaxis, decimals=0))
    ax.yaxis.set_major_formatter(comma_only_formatter(ax.yaxis, decimals=2))
    plt.tight_layout()
    plt.show()

    # plt.figure(figsize=(8, 4), facecolor='w', edgecolor='k')
    # plt.subplot(1, 2, 1)

    plt.figure(figsize=(8, 5), facecolor='w', edgecolor='k')
    gnlse.plot_wavelength_vs_distance_logarithmic(solution_, WL_range=[1000, 1100], cmap = _cmap_2d)

    # plt.subplot(1, 2, 2)
    # gnlse.plot_delay_vs_distance_logarithmic(solution_, time_range=[min(solution_.t),
    #                                                                 max(solution_.t)], cmap = _cmap_2d)

    ax = plt.gca()
    ax.xaxis.set_major_formatter(comma_only_formatter(ax.xaxis, decimals=0))
    ax.yaxis.set_major_formatter(comma_only_formatter(ax.yaxis, decimals=2))
    plt.tight_layout()
    plt.show()

    x_margin = 2
    plot_norm = None

    # _, ax = plt.subplots()
    # gnlse.plot_wavelength_for_distance_slice_logarithmic(solution_, WL_range=[3400, 3700], ax = ax, norm = plot_norm,
    #                                                      cmap="Blues", plot_for_zero=False)
    # ax.set_ylim([-190, -90])
    # ax.set_xlim([3400+x_margin, 3700-x_margin])
    # ax.grid(True)

    fig = plt.figure(figsize=(8, 5), constrained_layout=True)
    ax = fig.add_subplot(111)
    gnlse.plot_wavelength_for_distance_slice_logarithmic(solution_, WL_range=[3400, 3700], ax = ax,
                                                         z_slice = [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4],
                                                         norm = plot_norm, cmap=_cmap_1d,
                                                         plot_for_zero=True)
    ax.set_ylim([-210, -15])
    ax.set_xlim([3400+x_margin, 3700-x_margin])
    ax.grid(True)
    ax.xaxis.set_major_formatter(comma_only_formatter(ax.xaxis, decimals=0))
    ax.yaxis.set_major_formatter(comma_only_formatter(ax.yaxis, decimals=0))

    fig = plt.figure(figsize=(8, 5), constrained_layout=True)
    ax = fig.add_subplot(111)
    gnlse.plot_wavelength_for_distance_slice_logarithmic(solution_, WL_range=[620, 640], ax = ax,
                                                         z_slice = [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4],
                                                         norm = plot_norm, cmap=_cmap_1d,
                                                         plot_for_zero=True)
    ax.set_ylim([-210, 0])
    ax.set_xlim([620+x_margin, 640-x_margin])
    ax.grid(True)
    ax.xaxis.set_major_formatter(comma_only_formatter(ax.xaxis, decimals=0))
    ax.yaxis.set_major_formatter(comma_only_formatter(ax.yaxis, decimals=0))

    fig = plt.figure(figsize=(8, 5), constrained_layout=True)
    ax = fig.add_subplot(111)
    gnlse.plot_wavelength_for_distance_slice_logarithmic(solution_, WL_range=[1050, 1080], ax = ax,
                                                         z_slice = [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4],
                                                         norm = plot_norm, cmap=_cmap_1d,
                                                         plot_for_zero=True)
    ax.set_ylim([-210, 10])
    ax.set_xlim([1050+x_margin, 1080-x_margin])
    ax.grid(True)
    ax.xaxis.set_major_formatter(comma_only_formatter(ax.xaxis, decimals=0))
    ax.yaxis.set_major_formatter(comma_only_formatter(ax.yaxis, decimals=0))

    plt.legend(loc="upper right")
    plt.show()

    # ----------------------------------------

    fig = plt.figure(figsize=(8, 5), constrained_layout=True)
    ax = fig.add_subplot(111)
    gnlse.plot_wavelength_for_distance_slice_logarithmic(solution_, WL_range=[3400, 3700], ax = ax,
                                                         z_slice = [0, 0.03, 0.06, 0.09, 0.12, 0.15],
                                                         norm = plot_norm, cmap=_cmap_1d,
                                                         plot_for_zero=True)
    ax.set_ylim([-210, -90])
    ax.set_xlim([3400+x_margin, 3700-x_margin])
    ax.grid(True)
    ax.xaxis.set_major_formatter(comma_only_formatter(ax.xaxis, decimals=0))
    ax.yaxis.set_major_formatter(comma_only_formatter(ax.yaxis, decimals=0))

    fig = plt.figure(figsize=(8, 5), constrained_layout=True)
    ax = fig.add_subplot(111)
    gnlse.plot_wavelength_for_distance_slice_logarithmic(solution_, WL_range=[620, 640], ax = ax,
                                                         z_slice = [0, 0.03, 0.06, 0.09, 0.12, 0.15],
                                                         norm = plot_norm, cmap=_cmap_1d,
                                                         plot_for_zero=True)
    ax.set_ylim([-210, -85])
    ax.set_xlim([620+x_margin, 640-x_margin])
    ax.grid(True)
    ax.xaxis.set_major_formatter(comma_only_formatter(ax.xaxis, decimals=0))
    ax.yaxis.set_major_formatter(comma_only_formatter(ax.yaxis, decimals=0))

    fig = plt.figure(figsize=(8, 5), constrained_layout=True)
    ax = fig.add_subplot(111)
    gnlse.plot_wavelength_for_distance_slice_logarithmic(solution_, WL_range=[1050, 1080], ax = ax,
                                                         z_slice = [0, 0.03, 0.06, 0.09, 0.12, 0.15],
                                                         norm = plot_norm, cmap=_cmap_1d,
                                                         plot_for_zero=True)
    ax.set_ylim([-210, 10])
    ax.set_xlim([1050+x_margin, 1080-x_margin])
    ax.grid(True)
    ax.xaxis.set_major_formatter(comma_only_formatter(ax.xaxis, decimals=0))
    ax.yaxis.set_major_formatter(comma_only_formatter(ax.yaxis, decimals=0))

    plt.legend(loc="upper right")
    plt.show()

    # ----------------------------------------

    fig = plt.figure(figsize=(8, 5), constrained_layout=True)
    ax = fig.add_subplot(111)
    gnlse.visualization.plot_wavelength_slice_vs_distance_logarithmic(solution_,
                                                                      wavelengths=[626, 3550],
                                                                      ax = ax,
                                                                      norm = plot_norm, cmap=_cmap_1d)

    ax.set_xlim([0, max(solution_.Z)])
    ax.set_ylim([-250, 10])
    ax.grid(True)
    ax.xaxis.set_major_formatter(comma_only_formatter(ax.xaxis, decimals=2))
    ax.yaxis.set_major_formatter(comma_only_formatter(ax.yaxis, decimals=0))
    plt.show()


if __name__ == "__main__":
    path = os.path.join("solutions",
        "far_detuned_fwm" + \
        "_resolution_14" + \
        "_time_window_10" + \
        "_fiber_length_0.4" + \
        "_samples_1000" +\
        ".json")

    solution = gnlse.Solution()
    solution.from_file(path)
    plot_solution(solution)
    plot_gain(solution)