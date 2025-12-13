import os
import numpy as np
import gnlse
import math
from DFWM.plot_solutions import plot_solution
from DFWM.far_detuned_fwm import define_setup, solve_gnlse
from DFWM.plot_gain import plot_gain
import matplotlib as mpl
import matplotlib.pyplot as plt
import time, datetime

mpl.rcParams['font.family'] = 'Times New Roman'
mpl.rcParams['font.serif'] = ['Times New Roman']
mpl.rcParams['axes.unicode_minus'] = False

def calculate_average(solutions_list, *, weights=None):
    return gnlse.Solution(t=np.average([sol.t for sol in solutions_list],                        weights=weights, axis=0),
                          W=np.average([sol.W for sol in solutions_list],                        weights=weights, axis=0),
                          w_0=np.average([sol.w_0 for sol in solutions_list],                    weights=weights, axis=0),
                          Z=np.average([sol.Z for sol in solutions_list],                        weights=weights, axis=0),
                          At=np.sqrt(np.average([np.abs(sol.At) ** 2 for sol in solutions_list], weights=weights, axis=0)),
                          AW=np.sqrt(np.average([np.abs(sol.AW) ** 2 for sol in solutions_list], weights=weights, axis=0)))


def average_solutions(resolution, time_window, z_saves, wavelength, fiber_length, rtol, atol, peak_power,
                      n2, input_data_filepath, neff_max, samples):

    memory_usage_history = list()

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

    solutions_per_averaging = 2
    # ^ 3 mogłoby być lepsze, bo jest bliżej e, ale wydaje mi się że wtedy floating point errory będą większe

    solution_orders_num = math.ceil(math.log(samples, solutions_per_averaging))
    solutions_dict = {x: [] for x in range(solution_orders_num + 1)}

    time_start = time.time()
    print(f"Time Start: {datetime.datetime.fromtimestamp(time_start)}")


    try:
        for i in range(samples):
            print(f"{setup.pulse_model.name} {i + 1}/{samples}")
            solution = solve_gnlse(setup, n2=n2, json_path=input_data_filepath, neff_max=neff_max)
            solutions_dict[0].append(solution)

            for order_index, solutions_list in solutions_dict.items():

                if len(solutions_list) >= solutions_per_averaging:

                    average_solution = calculate_average(solutions_list)

                    solutions_dict[order_index] = []
                    solutions_dict[order_index+1].append(average_solution)

            solutions_mem_num = sum({key: len(value) for key, value in solutions_dict.items()}.values())
            print(f"Solutions in memory: {solutions_mem_num}")
            memory_usage_history.append((i+1, solutions_mem_num))
            time_now = time.time()
            print(f"Estimated finish time: {datetime.datetime.fromtimestamp((time_now - time_start) * (samples / (i+1)) + time_start)}")
    except KeyboardInterrupt:
        print("Process interrupted! (KeyboardInterrupt)")
    except MemoryError:
        print("Process interrupted! (MemoryError)")
    except:
        print("Process interrupted! (Unknown Error")

    solutions_with_weights =\
        [(calculate_average(value) if len(value) != 0 else None, len(value) * (solutions_per_averaging ** key))
         for key, value in solutions_dict.items()]

    solutions_with_weights = list(filter(lambda x: not(x[0] is None), solutions_with_weights))

    solutions, weights = zip(*solutions_with_weights)
    solution_averaged = calculate_average(solutions, weights=weights)

    path = os.path.join("solutions",
                        "far_detuned_fwm" + \
                        "_resolution_" + str(int(np.log(setup.resolution) / np.log(2))) + \
                        "_time_window_" + str(setup.time_window) + \
                        "_fiber_length_" + str(setup.fiber_length) + \
                        "_samples_" + str(samples) + \
                        ".json")

    solution_averaged.to_file(path)
    print(path)

    plot_solution(solution_averaged)
    plot_gain(solution_averaged)

    x, y = zip(*memory_usage_history)

    plt.plot(x, y)
    plt.xlabel("Liczba przeprowadzonych modelowań")
    plt.ylabel("Liczba modelowań w pamięci RAM")
    plt.grid(True)

    plt.show()

if __name__ == '__main__':

    average_solutions(resolution=2**14, time_window=10, z_saves=501, wavelength=1064, fiber_length=0.4, rtol=1e-3,
                      atol=1e-4, peak_power = 51000, n2=2.6e-20, neff_max=10, samples=1000,
                      input_data_filepath="../data/neff_far_detuned_FWM.json")
