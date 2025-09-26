import os
import numpy as np
import gnlse
from DFWM.plot_solutions import plot_solution


def define_setup(resolution, time_window, z_saves, wavelength, fiber_length,
                 raman_model, envelope, rtol, atol):
    setup_ = gnlse.GNLSESetup()

    # Numerical parameters
    setup_.resolution = resolution
    setup_.time_window = time_window  # ps
    setup_.z_saves = z_saves

    # Physical parameters
    setup_.wavelength = wavelength      # nm
    setup_.fiber_length = fiber_length  # m
    setup_.raman_model = raman_model    # gnlse.raman_blowwood
    setup_.self_steepening = True

    # Input pulse parameters
    peak_power = 51000  # W

    match envelope.lower():
        case "cw":       setup_.pulse_model = gnlse.CWEnvelope(peak_power, Pn=1e-16)
        case "gaussian": setup_.pulse_model = gnlse.GaussianEnvelope(peak_power, 1)
        case "sech":     setup_.pulse_model = gnlse.SechEnvelope(peak_power, 4)
        case _: raise ValueError

    if rtol is not None: setup_.rtol = 1e-8
    if atol is not None: setup_.atol = 1e-9

    return setup_

def solve_gnlse(setup_: gnlse.GNLSESetup):

    n2 = 2.6e-20  # m^2/W

    # read json file for neffs to cover interpolation example
    json_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'neff_far_detuned_FWM.json')
    json = gnlse.read_json(json_path)

    lambdas = json['neff'][:, 0] * 1e9 # wavelengths in nm
    neff = json['neff'][:, 1]          # efective refractive indices
    aeff = json['neff'][:, 2] * 1e-12  # efective mode area in m^2

    # The dispersion model is built from a Taylor expansion with coefficients given below.

    setup_.dispersion_model = gnlse.DispersionFiberFromInterpolation(loss=0.2, neff=neff, lambdas=lambdas,
                                                                    central_wavelength=setup_.wavelength)
    setup_.nonlinearity = gnlse.NonlinearityFromEffectiveArea(neff, aeff, lambdas, setup_.wavelength, n2=n2,
                                                              neff_max=10)

    print('%s...' % setup_.pulse_model.name)

    solver = gnlse.GNLSE(setup_)
    return solver.run()


if __name__ == '__main__':
    setup = define_setup(resolution=2**14, time_window=10, z_saves=51, wavelength=1064, fiber_length=0.2,
                         raman_model=None, envelope="CW", rtol=None, atol=None)

    solution = solve_gnlse(setup)
    plot_solution(solution)

    path = os.path.join("solutions",
           "far_detuned_fwm" + \
           "_resolution_" + str(int(np.log(setup.resolution) / np.log(2))) + \
           "_time_window_" + str(setup.time_window) + \
           "_fiber_length_" + str(setup.fiber_length) + \
           ".json")
    print(path)
    solution.to_file(path)
