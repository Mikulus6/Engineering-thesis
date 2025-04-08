import os
import numpy as np
import gnlse
from DFWM.far_detuned_fwm_plots import plot_solution

if __name__ == '__main__':
    setup = gnlse.GNLSESetup()

    setup.resolution = 2**14
    setup.time_window = 10  # ps
    setup.z_saves = 51

    setup.wavelength = 1064  # nm
    setup.fiber_length = 0.3  # m
    setup.raman_model = None #gnlse.raman_blowwood
    setup.self_steepening = True

    n2 = 2.6e-20  # m^2/W

    json_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'neff_far_detuned_FWM.json')
    json = gnlse.read_json(json_path)

    lambdas = json['neff'][:, 0] * 1e9
    neff = json['neff'][:, 1]
    Aeff = json['neff'][:, 2] * 1e-12

    loss = 0.2

    setup.dispersion_model = gnlse.DispersionFiberFromInterpolation(loss, neff, lambdas, setup.wavelength)
    setup.nonlinearity = gnlse.NonlinearityFromEffectiveArea(neff, Aeff, lambdas, setup.wavelength, n2=n2, neff_max=10)

    peak_power = 51000  # W

    samples = 200
    solutions = []

    for i in range(samples):
        setup.pulse_model = gnlse.CWEnvelope(peak_power, Pn=1e-16)
        print(f"{i+1}/{samples} {setup.pulse_model.name}...")
        solver = gnlse.GNLSE(setup)
        solution = solver.run()
        solutions.append(solution)

    average_solution = gnlse.Solution(t=np.mean([sol.t for sol in solutions], axis=0),
                                      W=np.mean([sol.W for sol in solutions], axis=0),
                                      w_0=np.mean([sol.w_0 for sol in solutions], axis=0),
                                      Z=np.mean([sol.Z for sol in solutions], axis=0),
                                      At=np.mean([np.abs(sol.At) for sol in solutions], axis=0),
                                      AW=np.mean([np.abs(sol.AW) for sol in solutions], axis=0))

    plot_solution(average_solution)

    path = "far_detuned_fwm" + \
           "_resolution_" + str(int(np.log(setup.resolution) / np.log(2))) + \
           "_time_window_" + str(setup.time_window) + \
           "_fiber_length_" + str(setup.fiber_length) + \
           "_samples_" + str(samples) + \
           ".json"

    print(path)
    average_solution.to_file(path)
