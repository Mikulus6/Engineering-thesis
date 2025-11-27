import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import gnlse

def plot_initial_parameters():

    mpl.rcParams['font.family'] = 'Times New Roman'
    mpl.rcParams['font.size'] = 12

    mpl.rcParams['figure.constrained_layout.use'] = True
    mpl.rcParams['figure.autolayout'] = True

    json_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'neff_far_detuned_FWM.json')
    json_data = gnlse.read_json(json_path)

    lambdas = json_data['neff'][:, 0] * 1e9
    neff = json_data['neff'][:, 1]
    aeff = json_data['neff'][:, 2] * 1e-2 # to plot easily (correct scalling is 1e-12)

    plt.figure()
    plt.plot(lambdas, neff, color="#0080FF")
    plt.xlabel("Długość fali [nm]")
    plt.ylabel("Efektywny współczynnik załamania")
    plt.grid(True)
    plt.tight_layout()

    plt.figure()
    plt.plot(lambdas, aeff, color="#0080FF")
    plt.xlabel("Długość fali [nm]")
    plt.ylabel("Efektywne pole modu [m² · 10⁻¹⁰]")
    plt.grid(True)
    plt.tight_layout()

    plt.show()

if __name__ == "__main__":
    plot_initial_parameters()