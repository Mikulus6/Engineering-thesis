import os
import numpy as np
import gnlse
from DFWM.plot_solutions import plot_solution
from DFWM.far_detuned_fwm import define_setup, solve_gnlse

if __name__ == '__main__':

    setup = define_setup(resolution=2**14, time_window=10, z_saves=51, wavelength=1064, fiber_length=0.4,
                         raman_model=None, envelope="CW", rtol=None, atol=None)

    samples = 200
    solutions = []

    for i in range(samples):
        print(f"{i+1}/{samples} {setup.pulse_model.name}...")
        solution = solve_gnlse(setup)
        solutions.append(solution)

    average_solution = gnlse.Solution(t=  np.mean([sol.t for sol in solutions],          axis=0),
                                      W=  np.mean([sol.W for sol in solutions],          axis=0),
                                      w_0=np.mean([sol.w_0 for sol in solutions],        axis=0),
                                      Z=  np.mean([sol.Z for sol in solutions],          axis=0),
                                      At= np.mean([np.abs(np.sqrt(sol.At)) for sol in solutions], axis=0)**2,
                                      AW= np.mean([np.abs(np.sqrt(sol.AW)) for sol in solutions], axis=0)**2)

    plot_solution(average_solution)

    path = os.path.join("solutions",
           "far_detuned_fwm" + \
           "_resolution_" + str(int(np.log(setup.resolution) / np.log(2))) + \
           "_time_window_" + str(setup.time_window) + \
           "_fiber_length_" + str(setup.fiber_length) + \
           "_samples_" + str(samples) + \
           ".json")

    print(path)
    average_solution.to_file(path)
