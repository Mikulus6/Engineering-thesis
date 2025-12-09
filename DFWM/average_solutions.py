import os
import numpy as np
import gnlse
from DFWM.plot_solutions import plot_solution
from DFWM.far_detuned_fwm import define_setup, solve_gnlse
from DFWM.plot_gain import plot_gain


def average_solutions(resolution, time_window, z_saves, wavelength, fiber_length, rtol, atol, peak_power,
                      n2, input_data_filepath, neff_max, samples):

    setup = define_setup(resolution=resolution,
                         time_window=time_window,
                         z_saves=z_saves,
                         wavelength=wavelength,
                         fiber_length=fiber_length,
                         raman_model=None,
                         envelope="CW",
                         rtol=rtol,
                         atol=atol,
                         peak_power=peak_power)
    solutions = []

    for i in range(samples):
        print(f"{i + 1}/{samples} {setup.pulse_model.name}...")
        solution = solve_gnlse(setup, n2=n2, json_path=input_data_filepath, neff_max=neff_max)
        solutions.append(solution)

    average_solution = gnlse.Solution(t=  np.mean([sol.t for sol in solutions],          axis=0),
                                      W=  np.mean([sol.W for sol in solutions],          axis=0),
                                      w_0=np.mean([sol.w_0 for sol in solutions],        axis=0),
                                      Z=  np.mean([sol.Z for sol in solutions],          axis=0),
                                      At= np.sqrt(np.mean([np.abs(sol.At)**2 for sol in solutions], axis=0)),
                                      AW= np.sqrt(np.mean([np.abs(sol.AW)**2 for sol in solutions], axis=0)))

    plot_solution(average_solution)
    plot_gain(average_solution)

    path = os.path.join("solutions",
                        "far_detuned_fwm" + \
                        "_resolution_" + str(int(np.log(setup.resolution) / np.log(2))) + \
                        "_time_window_" + str(setup.time_window) + \
                        "_fiber_length_" + str(setup.fiber_length) + \
                        "_samples_" + str(samples) + \
                        ".json")

    print(path)
    average_solution.to_file(path)

if __name__ == '__main__':

    average_solutions(resolution=2**14, time_window=10, z_saves=51, wavelength=1064, fiber_length=0.15, rtol=1e-3,
                      atol=1e-4, peak_power = 51000, n2=2.6e-20, neff_max=10, samples=1000,
                      input_data_filepath="../data/neff_far_detuned_FWM.json")
