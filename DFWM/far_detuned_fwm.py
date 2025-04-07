
import os

import numpy as np
import matplotlib.pyplot as plt

import gnlse


if __name__ == '__main__':
    setup = gnlse.GNLSESetup()

    # Numerical parameters
    setup.resolution = 2**14
    setup.time_window = 20  # ps

    setup.z_saves = 51

    # Physical parameters
    setup.wavelength = 10004  # nm
    setup.fiber_length = 0.1  # m
    setup.raman_model = None #gnlse.raman_blowwood
    setup.self_steepening = True

    #setup.rtol = 1e-8
    #setup.atol = 1e-9

    n2 = 2.6e-20  # m^2/W

    # read json file for neffs to cover interpolation example
    json_path = os.path.join(os.path.dirname(__file__), '..',
                             'data', 'neff_far_detuned_FWM.json')
    json = gnlse.read_json(json_path)

    # wavelengths in nm
    lambdas = json['neff'][:, 0] * 1e9
    # neffs
    neff = json['neff'][:, 1]
    # efective mode area in m^2
    Aeff = json['neff'][:, 2] * 1e-12

    # The dispersion model is built from a Taylor expansion with coefficients
    # given below.
    loss = 0.2

    setup.dispersion_model = gnlse.DispersionFiberFromInterpolation(
                                loss, neff, lambdas, setup.wavelength)

    setup.nonlinearity = gnlse.NonlinearityFromEffectiveArea(
        neff, Aeff, lambdas, setup.wavelength,
        n2=n2, neff_max=10)

    # Input pulse parameters
    peak_power = 51000  # W

    # This example extends the original code with additional simulations for
    setup.pulse_model = gnlse.CWEnvelope(peak_power, Pn=1e-16)
    # setup.pulse_model = gnlse.GaussianEnvelope(peak_power, 0.1)
    # setup.pulse_model = gnlse.SechEnvelope(peak_power, 4)

    print('%s...' % setup.pulse_model.name)

    solver = gnlse.GNLSE(setup)
    solution = solver.run()

    plt.figure(figsize=(8, 4), facecolor='w', edgecolor='k')
    plt.subplot(1, 2, 1)
    gnlse.plot_wavelength_vs_distance(
        solution,
        WL_range=[300, 4000],
        cmap = 'jet')

    plt.subplot(1, 2, 2)
    gnlse.plot_delay_vs_distance(
        solution,
        time_range=[-10, 10],
        cmap = 'jet')

    plt.tight_layout()
    plt.show()

    plt.figure(figsize=(8, 4), facecolor='w', edgecolor='k')
    plt.subplot(1, 2, 1)
    gnlse.plot_wavelength_vs_distance(
        solution,
        WL_range=[1000, 1100],
        cmap = 'jet')

    plt.subplot(1, 2, 2)
    gnlse.plot_delay_vs_distance(
        solution,
        time_range=[-5, 5],
        cmap = 'jet')

    plt.tight_layout()
    plt.show()


    plt.figure(2, figsize=(8,4), dpi=300)
    fig2, ax = plt.subplots()
    gnlse.plot_wavelength_for_distance_slice_logarithmic(
        solution,
        WL_range=[3400, 3700],
        ax = ax,
        norm = 51e3)
    ax.set_ylim([-100, 5])


    plt.figure(2, figsize=(8,4), dpi=300)
    fig2, ax = plt.subplots()
    gnlse.plot_wavelength_for_distance_slice_logarithmic(
        solution,
        WL_range=[3400, 3700],
        ax = ax,
        z_slice = [0., 0.01, 0.05, .10],
        norm = 51e3)
    ax.set_ylim([-400, 5])
    ax.grid(True)

    plt.figure(2, figsize=(8,4), dpi=300)
    fig2, ax = plt.subplots()
    gnlse.plot_wavelength_for_distance_slice_logarithmic(
        solution,
        WL_range=[620, 640],
        ax = ax,
        z_slice = [0., 0.01, 0.05, .10],
        norm = 51e3)
    ax.set_ylim([-400, 5])
    ax.grid(True)

    path = "far_detuned_fwm" + \
        "_resolution_" + str(int(np.log(setup.resolution)/np.log(2))) + \
        "_time_window_" + str(setup.time_window) + \
        "_fiber_length_" + str(setup.fiber_length) + \
        ".json"
    print(path)
    solution.to_file(path)
