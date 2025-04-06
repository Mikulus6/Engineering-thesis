import matplotlib.pyplot as plt
import gnlse


if __name__ == "__main__":
    path = "far_detuned_fwm" + \
        "_resolution_14" + \
        "_time_window_20" + \
        "_fiber_length_0.5" + \
        ".json"

    solution = gnlse.Solution()
    solution.from_file(path)

    plt.figure(figsize=(8, 4), facecolor='w', edgecolor='k')
    plt.subplot(1, 2, 1)
    gnlse.plot_wavelength_vs_distance_logarithmic(
        solution, 
        WL_range=[300, 4000],
        cmap = 'jet')

    plt.subplot(1, 2, 2)
    gnlse.plot_delay_vs_distance_logarithmic(
        solution, 
        time_range=[-40, 40],
        cmap = 'jet')

    plt.tight_layout()
    plt.show()
    
    plt.figure(figsize=(8, 4), facecolor='w', edgecolor='k')
    plt.subplot(1, 2, 1)
    gnlse.plot_wavelength_vs_distance_logarithmic(
        solution, 
        WL_range=[1000, 1100],
        cmap = 'jet')

    plt.subplot(1, 2, 2)
    gnlse.plot_delay_vs_distance_logarithmic(
        solution, 
        time_range=[-40, 40],
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
        z_slice = solution.Z[0::10],
        norm = 51e3)
    ax.set_ylim([-200, 50])
    ax.grid(True)
    
    plt.figure(2, figsize=(8,4), dpi=300)
    fig2, ax = plt.subplots()
    gnlse.plot_wavelength_for_distance_slice_logarithmic(
        solution, 
        WL_range=[620, 640],
        ax = ax,
        z_slice = solution.Z[0::10],
        norm = 51e3)
    ax.set_ylim([-200, 50])
    ax.grid(True)
    
    plt.figure(2, figsize=(8,4), dpi=300)
    fig2, ax = plt.subplots()
    gnlse.plot_wavelength_for_distance_slice_logarithmic(
        solution, 
        WL_range=[1050, 1080],
        ax = ax,
        z_slice = solution.Z[0::10],
        norm = 51e3)
    ax.set_ylim([-200, 50])
    ax.grid(True)
