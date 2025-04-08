import os
import numpy as np
import gnlse
from DFWM.far_detuned_fwm_plots import plot_solution

if __name__ == '__main__':
    setup = gnlse.GNLSESetup()

    # Numerical parameters
    setup.resolution = 2**14
    setup.time_window = 10  # ps

    setup.z_saves = 51

    # Physical parameters
    setup.wavelength = 1064  # nm
    setup.fiber_length = 0.2  # m
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
    # setup.pulse_model = gnlse.GaussianEnvelope(peak_power, 1)
    # setup.pulse_model = gnlse.SechEnvelope(peak_power, 4)

    print('%s...' % setup.pulse_model.name)

    solver = gnlse.GNLSE(setup)
    solution = solver.run()

    plot_solution(solution)

    path = "far_detuned_fwm" + \
        "_resolution_" + str(int(np.log(setup.resolution)/np.log(2))) + \
        "_time_window_" + str(setup.time_window) + \
        "_fiber_length_" + str(setup.fiber_length) + \
        ".json"
    print(path)
    solution.to_file(path)
