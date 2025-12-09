"""
Various plotting functions for visualizing GNLSE simulations using Matplotlib.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import TwoSlopeNorm, Normalize

from scipy.interpolate import RectBivariateSpline

from gnlse.common import c


def plot_frequency_vs_distance_logarithmic(solver, ax=None, norm=None,
                                           frequency_range=None, cmap="magma"):
    """Plotting results in logarithmic scale in frequency domain.

    Parameters
    ----------
    solver : Solution
        Model outputs in the form of a ``Solution`` object.
    frequency_range : list, (2, )
        frequency range. Set [-150, 150] as default.
    ax : :class:`~matplotlib.axes.Axes`
        :class:`~matplotlib.axes.Axes` instance for plotting
    norm : float
        Normalization factor for output spectrum. As default maximum of
        square absolute of ``solver.AW`` variable is taken.

    Returns
    -------
    ax : :class:`~matplotlib.axes.Axes`
       Used :class:`~matplotlib.axes.Axes` instance.
    """

    if ax is None:
        ax = plt.gca()

    if norm is None:
        norm = np.max(np.abs(solver.AW)**2)

    lIW = np.fliplr(
        10 * np.log10(np.abs(solver.AW)**2 / norm,
                      where=(np.abs(solver.AW)**2 > 0)))
    frequency = (solver.W - solver.w_0) / 2 / np.pi  # frequency grid

    if frequency_range is not None:
        iis = np.logical_and(frequency >= frequency_range[0],
                             frequency <= frequency_range[1])
        # indices of interest

        frequency = frequency[iis]
        lIW = lIW[:, iis]

    ax.imshow(lIW, origin='lower', aspect='auto', cmap=cmap,
              extent=[np.min(frequency), np.max(frequency),
                      0, np.max(solver.Z)], vmin=-40)
    ax.set_xlabel("Frequency [THz]")
    ax.set_ylabel("Distance [m]")
    return ax


def plot_frequency_vs_distance(solver, frequency_range=None,
                               ax=None, norm=None, cmap="magma"):
    """Plotting results in frequency domain. Linear scale.

    Parameters
    ----------
    solver : Solution
        Model outputs in the form of a ``Solution`` object.
    frequency_range : list, (2, )
        frequency range. Set [-150, 150] as default.
    ax : :class:`~matplotlib.axes.Axes`
        :class:`~matplotlib.axes.Axes` instance for plotting
    norm : float
        Normalization factor for output spectrum. As default maximum of
        square absolute of ``solver.AW`` variable is taken.

    Returns
    -------
    ax : :class:`~matplotlib.axes.Axes`
       Used :class:`~matplotlib.axes.Axes` instance.
    """

    if ax is None:
        ax = plt.gca()

    if norm is None:
        norm = np.max(np.abs(solver.AW)**2)

    IW = np.fliplr(
        np.abs(solver.AW)**2 / norm)
    frequency = (solver.W - solver.w_0) / 2 / np.pi  # frequency grid

    if frequency_range is not None:
        iis = np.logical_and(frequency >= frequency_range[0],
                             frequency <= frequency_range[1])
        # indices of interest

        frequency = frequency[iis]
        IW = IW[:, iis]

    ax.imshow(IW, origin='lower', aspect='auto', cmap=cmap,
              extent=[np.min(frequency), np.max(frequency),
                      0, np.max(solver.Z)], vmin=0)
    ax.set_xlabel("Frequency [THz]")
    ax.set_ylabel("Distance [m]")
    return ax


def plot_delay_for_distance_slice(solver, time_range=None, ax=None,
                                  z_slice=None, norm=None):
    """Plotting intensity in linear scale in time domain.

    Parameters
    ----------
    solver : Solution
        Model outputs in the form of a ``Solution`` object.
    time_range : list, (2, )
        Time range. Set [min(``solver.t``), max(``solver.t``)] as default.
    ax : :class:`~matplotlib.axes.Axes`
        :class:`~matplotlib.axes.Axes` instance for plotting.
    norm : float
        Normalization factor for output spectrum. As default maximum of
        square absolute of ``solver.At`` variable is taken.

    Returns
    -------
    ax : :class:`~matplotlib.axes.Axes`
       Used :class:`~matplotlib.axes.Axes` instance.
    """

    if ax is None:
        ax = plt.gca()

    if time_range is None:
        time_range = [np.min(solver.t), np.max(solver.t)]

    if norm is None:
        norm = np.max(np.abs(solver.At[0][:])**2)

    It = np.abs(solver.At)**2 / norm

    # indices of interest if no z_slice positions were given
    if z_slice is None:
        iis = [0, -1]
    # indices of interest nearest to given z_slice positions
    else:
        iis = [np.nonzero(
            np.min(np.abs(solver.Z - z)) == np.abs(solver.Z - z)
        )[0][0] for z in z_slice]

    for i in iis:
        label_i = "z = " + str(solver.Z[i]) + "m"
        ax.plot(solver.t, It[i][:], label=label_i)

    ax.set_xlim(time_range)
    ax.set_xlabel("Delay [ps]")
    ax.set_ylabel("Normalized Power")
    ax.legend()
    return ax


def plot_wavelength_for_distance_slice(solver, WL_range=None, ax=None,
                                       z_slice=None, norm=None):
    """Plotting chosen slices of intensity
    in linear scale in wavelength domain.

    Parameters
    ----------
    solver : Solution
        Model outputs in the form of a ``Solution`` object.
    WL_range : list, (2, )
        Wavelength range. Set [400, 1350] as default.
    ax : :class:`~matplotlib.axes.Axes`
        :class:`~matplotlib.axes.Axes` instance for plotting
    norm : float
        Normalization factor for output spectrum. As default maximum of
        square absolute of ``solver.AW`` variable is taken.

    Returns
    -------
    ax : :class:`~matplotlib.axes.Axes`
       Used :class:`~matplotlib.axes.Axes` instance.
    """

    if ax is None:
        ax = plt.gca()

    if WL_range is None:
        WL_range = [np.min(c / (solver.W / 2 / np.pi)),
                    np.max(c / (solver.W / 2 / np.pi))]

    if norm is None:
        norm = np.max(np.abs(solver.AW)**2)

    IW = np.fliplr(
        np.abs(solver.AW)**2 / norm)
    WL = 2 * np.pi * c / solver.W  # wavelength grid
    WL_asc = np.flip(WL, )  # ascending order for interpolation
    iio = np.logical_and(WL_asc > WL_range[0],
                         WL_asc < WL_range[1])  # indices in order

    WL_asc = WL_asc[iio]
    IW = IW[:, iio]

    # indices of interest if no z_slice positions were given
    if z_slice is None:
        iis = [0, -1]
    # indices of interest nearest to given z_slice positions
    else:
        iis = [np.nonzero(
            np.min(np.abs(solver.Z - z)) == np.abs(solver.Z - z)
        )[0][0] for z in z_slice]

    for i in iis:
        label_i = "z = " + str(solver.Z[i]) + "m"
        ax.plot(WL_asc, IW[i][:], label=label_i)

    ax.set_xlim([np.min(WL_asc), np.max(WL_asc)])
    ax.set_xlabel("Wavelength [nm]")
    ax.set_ylabel("Normalized Spectral Density")
    ax.legend()
    return ax


def plot_wavelength_for_distance_slice_logarithmic(solver, WL_range=None,
                                                   ax=None,
                                                   z_slice=None, norm=None, cmap="magma",
                                                   plot_gain=False,
                                                   plot_for_zero=False):
    """Plotting chosen slices of intensity
    in linear scale in wavelength domain.

    Parameters
    ----------
    solver : Solution
        Model outputs in the form of a ``Solution`` object.
    WL_range : list, (2, )
        Wavelength range. Set [400, 1350] as default.
    ax : :class:`~matplotlib.axes.Axes`
        :class:`~matplotlib.axes.Axes` instance for plotting
    norm : float
        Normalization factor for output spectrum. As default maximum of
        square absolute of ``solver.AW`` variable is taken.

    Returns
    -------
    ax : :class:`~matplotlib.axes.Axes`
       Used :class:`~matplotlib.axes.Axes` instance.
    """

    if ax is None:
        ax = plt.gca()

    if WL_range is None:
        WL_range = [np.min(c / (solver.W / 2 / np.pi)),
                    np.max(c / (solver.W / 2 / np.pi))]

    if norm is None:
        norm = np.max(np.abs(solver.AW)**2)

    lIW = np.fliplr(
        10 * np.log10(np.abs(solver.AW)**2 / norm,
                      where=(np.abs(solver.AW)**2 > 0)))
    WL = 2 * np.pi * c / solver.W  # wavelength grid
    WL_asc = np.flip(WL, )  # ascending order for interpolation
    iio = np.logical_and(WL_asc > WL_range[0],
                         WL_asc < WL_range[1])  # indices in order

    WL_asc = WL_asc[iio]
    lIW = lIW[:, iio]

    if plot_gain:
        lIW2 = np.zeros_like(lIW)
        for x in range(lIW.shape[1]):
            for y in range(lIW.shape[0]):
                try:
                    lIW2[y, x] = (lIW[y+1, x] - lIW[y, x])/(max(solver.Z)/len(solver.Z))
                except IndexError:
                    lIW2[y, x] = lIW2[y-1, x]

        lIW = lIW2

    # indices of interest if no z_slice positions were given
    if z_slice is None:
        iis = [0, -1]
    # indices of interest nearest to given z_slice positions
    else:
        iis = [np.nonzero(
            np.min(np.abs(solver.Z - z)) == np.abs(solver.Z - z)
        )[0][0] for z in z_slice]

    cmap_obj = cm.get_cmap(cmap)
    for index_, i in enumerate(iis):

        if not plot_for_zero and index_ == 0:
            continue

        label_i = "z = " + str(round(solver.Z[i], ndigits=3)) + "m"
        ax.plot(WL_asc, lIW[i][:], label=label_i, color=cmap_obj(index_/(len(iis)-1)))

    ax.set_ylim(-40)
    ax.set_xlabel("Długość fali [nm]")
    ax.set_ylabel("Znormalizowana gęstość spektralna" if not plot_gain else "Przyrost gęstości spektralnej [dB/m]")
    ax.legend()
    return ax


def plot_delay_for_distance_slice_logarithmic(solver, time_range=None, ax=None,
                                              z_slice=None, norm=None):
    """Plotting chosen slices of intensity in linear scale in time domain.

    Parameters
    ----------
    solver : Solution
        Model outputs in the form of a ``Solution`` object.
    time_range : list, (2, )
        Time range. Set [min(``solver.t``), max(``solver.t``)] as default.
    ax : :class:`~matplotlib.axes.Axes`
        :class:`~matplotlib.axes.Axes` instance for plotting.
    norm : float
        Normalization factor for output spectrum. As default maximum of
        square absolute of ``solver.At`` variable is taken.

    Returns
    -------
    ax : :class:`~matplotlib.axes.Axes`
       Used :class:`~matplotlib.axes.Axes` instance.
    """

    if ax is None:
        ax = plt.gca()

    if time_range is None:
        time_range = [np.min(solver.t), np.max(solver.t)]

    if norm is None:
        norm = np.max(np.abs(solver.At)**2)

    lIt = 10 * np.log10(np.abs(solver.At)**2 / norm,
                        where=(np.abs(solver.At)**2 > 0))

    # indices of interest if no z_slice positions were given
    if z_slice is None:
        iis = [0, -1]
    # indices of interest nearest to given z_slice positions
    else:
        iis = [np.nonzero(
            np.min(np.abs(solver.Z - z)) == np.abs(solver.Z - z)
        )[0][0] for z in z_slice]

    for i in iis:
        label_i = "z = " + str(solver.Z[i]) + "m"
        ax.plot(solver.t, lIt[i][:], label=label_i)

    ax.set_xlim(time_range)
    ax.set_ylim(-40)
    ax.set_xlabel("Delay [ps]")
    ax.set_ylabel("Normalized Power")
    ax.legend()
    return ax


def plot_frequency_for_distance_slice(solver, frequency_range=None, ax=None,
                                      z_slice=None, norm=None):
    """Plotting chosen slices of intensity in linear scale in frequency domain.

    Parameters
    ----------
    solver : Solution
        Model outputs in the form of a ``Solution`` object.
    frequency_range : list, (2, )
        frequency range. Set [-150, 150] as default.
    ax : :class:`~matplotlib.axes.Axes`
        :class:`~matplotlib.axes.Axes` instance for plotting
    norm : float
        Normalization factor for output spectrum. As default maximum of
        square absolute of ``solver.AW`` variable is taken.

    Returns
    -------
    ax : :class:`~matplotlib.axes.Axes`
       Used :class:`~matplotlib.axes.Axes` instance.
    """

    if ax is None:
        ax = plt.gca()

    if frequency_range is None:
        frequency_range = [np.min((solver.W - solver.w_0) / 2 / np.pi),
                           np.max((solver.W - solver.w_0) / 2 / np.pi)]

    if norm is None:
        norm = np.max(np.abs(solver.AW)**2)

    IW = np.fliplr(
        np.abs(solver.AW)**2 / norm)

    # indices of interest if no z_slice positions were given
    if z_slice is None:
        iis = [0, -1]  # beginning, end
    # indices of interest nearest to given z_slice positions
    else:
        iis = [np.nonzero(
            np.min(np.abs(solver.Z - z)) == np.abs(solver.Z - z)
        )[0][0] for z in z_slice]

    for i in iis:
        label_i = "z = " + str(solver.Z[i]) + "m"
        ax.plot((solver.W - solver.w_0) / 2 / np.pi, IW[i][:], label=label_i)

    ax.set_xlim(frequency_range)
    ax.set_xlabel("Frequency [Thz]")
    ax.set_ylabel("Normalized Spectral Density")
    ax.legend()
    return ax


def plot_frequency_for_distance_slice_logarithmic(solver, frequency_range=None,
                                                  ax=None, z_slice=None,
                                                  norm=None):
    """Plotting chosen slices of intensity
    in logarithmic scale in frequency domain.

    Parameters
    ----------
    solver : Solution
        Model outputs in the form of a ``Solution`` object.
    frequency_range : list, (2, )
        frequency range. Set [-150, 150] as default.
    ax : :class:`~matplotlib.axes.Axes`
        :class:`~matplotlib.axes.Axes` instance for plotting
    norm : float
        Normalization factor for output spectrum. As default maximum of
        square absolute of ``solver.AW`` variable is taken.

    Returns
    -------
    ax : :class:`~matplotlib.axes.Axes`
       Used :class:`~matplotlib.axes.Axes` instance.
    """

    if ax is None:
        ax = plt.gca()

    if frequency_range is None:
        frequency_range = [np.min((solver.W - solver.w_0) / 2 / np.pi),
                           np.max((solver.W - solver.w_0) / 2 / np.pi)]

    if norm is None:
        norm = np.max(np.abs(solver.AW)**2)

    lIW = np.fliplr(
        10 * np.log10(np.abs(solver.AW)**2 / norm,
                      where=(np.abs(solver.AW)**2 > 0)))

    # indices of interest if no z_slice positions were given
    if z_slice is None:
        iis = [0, -1]  # beginning, end

    # indices of interest nearest to given z_slice positions
    else:
        iis = [np.nonzero(
            np.min(np.abs(solver.Z - z)) == np.abs(solver.Z - z)
        )[0][0] for z in z_slice]

    for i in iis:
        label_i = "z = " + str(solver.Z[i]) + "m"
        ax.plot((solver.W - solver.w_0) / 2 / np.pi, lIW[i][:], label=label_i)

    ax.set_xlim(frequency_range)
    ax.set_ylim(-40)
    ax.set_xlabel("Frequency [Thz]")
    ax.set_ylabel("Normalized Spectral Density")
    ax.legend()
    return ax


def plot_delay_vs_distance_logarithmic(solver, time_range=None, ax=None,
                                       norm=None, cmap="magma"):
    """Plotting intensity in logarithmic scale in time domain.

    Parameters
    ----------
    solver : Solution
        Model outputs in the form of a ``Solution`` object.
    time_range : list, (2, )
        Time range. Set [min(``solver.t``), max(``solver.t``)] as default.
    ax : :class:`~matplotlib.axes.Axes`
        :class:`~matplotlib.axes.Axes` instance for plotting.
    norm : float
        Normalization factor for output spectrum. As default maximum of
        square absolute of ``solver.At`` variable is taken.

    Returns
    -------
    ax : :class:`~matplotlib.axes.Axes`
       Used :class:`~matplotlib.axes.Axes` instance.
    """
    if ax is None:
        ax = plt.gca()

    if time_range is None:
        time_range = [np.min(solver.t), np.max(solver.t)]

    if norm is None:
        norm = np.max(np.abs(solver.At)**2)

    lIT = 10 * np.log10(np.abs(solver.At)**2 / norm,
                        where=(np.abs(solver.At)**2 > 0))

    ax.pcolormesh(solver.t, solver.Z, lIT, shading="auto", vmin=-40,
                  cmap=cmap)
    ax.set_xlim(time_range)
    ax.set_xlabel("Delay [ps]")
    ax.set_ylabel("Distance [m]")
    return ax


def plot_delay_vs_distance(solver, time_range=None, ax=None, norm=None,
                           cmap="magma"):
    """Plotting normalized intensity in linear scale in time domain.

    Parameters
    ----------
    solver : Solution
        Model outputs in the form of a ``Solution`` object.
    time_range : list, (2, )
        Time range. Set [min(``solver.t``), max(``solver.t``)] as default.
    ax : :class:`~matplotlib.axes.Axes`
        :class:`~matplotlib.axes.Axes` instance for plotting.
    norm : float
        Normalization factor for output spectrum. As default maximum of
        square absolute of ``solver.At`` variable is taken.

    Returns
    -------
    ax : :class:`~matplotlib.axes.Axes`
       Used :class:`~matplotlib.axes.Axes` instance.
    """
    if ax is None:
        ax = plt.gca()

    if time_range is None:
        time_range = [np.min(solver.t), np.max(solver.t)]

    if norm is None:
        norm = np.max(np.abs(solver.At)**2)

    lIT = np.abs(solver.At)**2 / norm

    ax.pcolormesh(solver.t, solver.Z, lIT, shading="auto", vmin=0,
                  cmap=cmap)
    ax.set_xlim(time_range)
    ax.set_xlabel("Delay [ps]")
    ax.set_ylabel("Distance [m]")
    return ax


def plot_wavelength_vs_distance(solver, WL_range=None, ax=None,
                                norm=None, cmap="magma"):
    """Plotting results in linear scale in wavelength domain.

    Parameters
    ----------
    solver : Solution
        Model outputs in the form of a ``Solution`` object.
    WL_range : list, (2, )
        Wavelength range. Set [400, 1350] as default.
    ax : :class:`~matplotlib.axes.Axes`
        :class:`~matplotlib.axes.Axes` instance for plotting
    norm : float
        Normalization factor for output spectrum. As default maximum of
        square absolute of ``solver.AW`` variable is taken.

    Returns
    -------
    ax : :class:`~matplotlib.axes.Axes`
       Used :class:`~matplotlib.axes.Axes` instance.
    """

    if ax is None:
        ax = plt.gca()

    if WL_range is None:
        WL_range = [np.min(c / (solver.W / 2 / np.pi)),
                    np.max(c / (solver.W / 2 / np.pi))]

    if norm is None:
        norm = np.max(np.abs(solver.AW)**2)

    IW = np.fliplr(
        np.abs(solver.AW)**2 / norm)
    WL = 2 * np.pi * c / solver.W  # wavelength grid
    WL_asc = np.flip(WL, )  # ascending order for interpolation
    iis = np.logical_and(WL_asc > WL_range[0],
                         WL_asc < WL_range[1])  # indices of interest

    WL_asc = WL_asc[iis]
    IW = IW[:, iis]

    interpolator = RectBivariateSpline(solver.Z, WL_asc, IW)
    newWL = np.linspace(np.min(WL_asc), np.max(WL_asc), IW.shape[1])
    toshow = interpolator(solver.Z, newWL)

    ax.imshow(toshow, origin='lower', aspect='auto', cmap=cmap,
              extent=[np.min(WL_asc), np.max(WL_asc), 0, np.max(solver.Z)])
    ax.set_xlabel("Wavelength [nm]")
    ax.set_ylabel("Distance [m]")
    return ax


def plot_wavelength_vs_distance_logarithmic(solver, WL_range=None,
                                            ax=None, norm=None, cmap="magma",
                                            *, plot_gain=False, use_zero_norm: bool = False,
                                            norm_scale_inverse_factor: int | float = 1):
    """Plotting results in logarithmic scale in wavelength domain.

    Parameters
    ----------
    solver : Solution
        Model outputs in the form of a ``Solution`` object.
    WL_range : list, (2, )
        Wavelength range. Set [400, 1350] as default.
    ax : :class:`~matplotlib.axes.Axes`
        :class:`~matplotlib.axes.Axes` instance for plotting
    norm : float
        Normalization factor for output spectrum. As default maximum of
        square absolute of ``solver.AW`` variable is taken.

    Returns
    -------
    ax : :class:`~matplotlib.axes.Axes`
       Used :class:`~matplotlib.axes.Axes` instance.
    """

    if ax is None:
        ax = plt.gca()

    if WL_range is None:
        WL_range = [np.min(c / (solver.W / 2 / np.pi)),
                    np.max(c / (solver.W / 2 / np.pi))]

    if norm is None:
        norm = np.max(np.abs(solver.AW)**2)

    lIW = np.fliplr(
        10 * np.log10(np.abs(solver.AW)**2 / norm,
                      where=(np.abs(solver.AW)**2 > 0)))
    WL = 2 * np.pi * c / solver.W  # wavelength grid
    WL_asc = np.flip(WL, )  # ascending order for interpolation
    iis = np.logical_and(WL_asc > WL_range[0],
                         WL_asc < WL_range[1])  # indices of interest

    WL_asc = WL_asc[iis]
    lIW = lIW[:, iis]

    if plot_gain:
        lIW2 = np.zeros_like(lIW)
        for x in range(lIW.shape[1]):
            for y in range(lIW.shape[0]):
                try:
                    lIW2[y, x] = (lIW[y+1, x] - lIW[y, x])/(max(solver.Z)/len(solver.Z))
                except IndexError:
                    lIW2[y, x] = lIW2[y-1, x]
            print(f"Derivating {x}/{lIW.shape[1]}")
        lIW = lIW2

    interpolator = RectBivariateSpline(solver.Z, WL_asc, lIW)
    newWL = np.linspace(np.min(WL_asc), np.max(WL_asc), lIW.shape[1])
    toshow = interpolator(solver.Z, newWL)

    ax.set_xlabel("Długość fali [nm]")
    ax.set_ylabel("Długość światłowodu [m]")

    if not use_zero_norm:
        norm_obj = Normalize(vmin=np.nanmin(lIW), vmax=np.nanmax(lIW))
    else:
        value = max(abs(np.nanmin(lIW)), abs(np.nanmax(lIW)))
        norm_obj = TwoSlopeNorm(vmin=-value/norm_scale_inverse_factor, vcenter=0, vmax=value/norm_scale_inverse_factor)

    ax.imshow(toshow, origin='lower', aspect='auto', cmap=cmap,
              norm = norm_obj,
              extent=[np.min(WL_asc), np.max(WL_asc), 0, np.max(solver.Z)])

    sm = cm.ScalarMappable(norm=norm_obj, cmap=cmap)
    sm.set_array([])

    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label("Znormalizowana gęstość spektralna" if not plot_gain else "Przyrost gęstości spektralnej [dB/m]")

    return ax


def plot_wavelength_slice_vs_distance_logarithmic(solver, wavelengths,
                                                 ax=None, norm=None, cmap="magma",
                                                 *, plot_gain=False):
    """

    Parameters
    ----------
    solver : Solution
        Model outputs in the form of a ``Solution`` object.
    wavelengths:
        Wavelengths list with values for slices
    ax : :class:`~matplotlib.axes.Axes`
        :class:`~matplotlib.axes.Axes` instance for plotting
    norm : float
        Normalization factor for output spectrum. As default maximum of
        square absolute of ``solver.AW`` variable is taken.

    Returns
    -------
    ax : :class:`~matplotlib.axes.Axes`
       Used :class:`~matplotlib.axes.Axes` instance.
    """

    if isinstance(wavelengths, int) or isinstance(wavelengths, float):
        wavelengths = [wavelengths]

    if ax is None:
        ax = plt.gca()

    WL_range = [min(wavelengths)-1, max(wavelengths)+1]

    if norm is None:
        norm = np.max(np.abs(solver.AW)**2)

    lIW = np.fliplr(
        10 * np.log10(np.abs(solver.AW)**2 / norm,
                      where=(np.abs(solver.AW)**2 > 0)))
    WL = 2 * np.pi * c / solver.W  # wavelength grid
    WL_asc = np.flip(WL, )  # ascending order for interpolation
    iis = np.logical_and(WL_asc > WL_range[0],
                         WL_asc < WL_range[1])  # indices of interest

    WL_asc = WL_asc[iis]
    lIW = lIW[:, iis]

    if plot_gain:
        lIW2 = np.zeros_like(lIW)
        for x in range(lIW.shape[1]):
            for y in range(lIW.shape[0]):
                try:
                    lIW2[y, x] = (lIW[y+1, x] - lIW[y, x])/(max(solver.Z)/len(solver.Z))
                except IndexError:
                    lIW2[y, x] = lIW2[y-1, x]

        lIW = lIW2

    interpolator = RectBivariateSpline(solver.Z, WL_asc, lIW)
    newWL = np.linspace(np.min(WL_asc), np.max(WL_asc), lIW.shape[1])
    toshow = interpolator(solver.Z, newWL)

    ax.set_xlabel("Długość światłowodu [m]")
    ax.set_ylabel("Znormalizowana gęstość spektralna" if not plot_gain else "Przyrost gęstości spektralnej [dB/m]")

    cmap_obj = cm.get_cmap(cmap)
    for index_, wavelength in enumerate(wavelengths):
        idx = np.argmin(np.abs(newWL - wavelength))

        # Extract slice vs propagation distance Z
        slice_vals = toshow[:, idx]

        # Plot the slice
        ax.plot(solver.Z, slice_vals, label=f"λ={wavelength}nm",
                color=cmap_obj(index_/max((len(wavelengths)-1), 1)))

    ax.legend()

    return ax


def quick_plot(solution):
    """Plotting results in time and frequency domain for default value
    of parameters.

    Parameters
    ----------
    solver : Solution
        Model outputs in the form of a ``Solution`` object.
    """
    plt.suptitle('GNLSE solution')

    plt.subplot(1, 2, 1)
    plot_wavelength_vs_distance(solution)

    plt.subplot(1, 2, 2)
    plot_delay_vs_distance(solution)

    plt.show()
