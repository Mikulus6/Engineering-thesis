import os
import matplotlib.pyplot as plt
import gnlse

def plot_initial_parameters():

    json_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'neff_far_detuned_FWM.json')
    json_data = gnlse.read_json(json_path)

    lambdas = json_data['neff'][:, 0] * 1e9  # wavelengths in nm
    neff = json_data['neff'][:, 1]           # effective refractive indices
    aeff = json_data['neff'][:, 2] * 1e-12   # effective mode area in m^2

    plt.figure()
    plt.plot(lambdas, neff)
    plt.xlabel("Wavelength (nm)")
    plt.ylabel("Effective Refractive Index (n_eff)")
    plt.grid(True)

    plt.figure()
    plt.plot(lambdas, aeff)
    plt.xlabel("Wavelength (nm)")
    plt.ylabel("Effective Mode Area (m²)")
    plt.grid(True)

    plt.show()

if __name__ == "__main__":
    plot_initial_parameters()