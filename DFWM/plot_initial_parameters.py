import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import gnlse
from matplotlib.ticker import FuncFormatter

def comma_only_formatter(axis, decimals=2):
    base = axis.get_major_formatter()
    def fmt(x, pos):
        s = base.format_data_short(x)
        return f"{float(s):.{decimals}f}".replace('.', ',') if decimals > 0 else str(int(round(float(s))))
    return FuncFormatter(fmt)

def plot_initial_parameters():

    mpl.rcParams['font.family'] = 'Times New Roman'
    mpl.rcParams['font.size'] = 16

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
    ax = plt.gca()
    ax.yaxis.set_major_formatter(comma_only_formatter(ax.yaxis, decimals=2))
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