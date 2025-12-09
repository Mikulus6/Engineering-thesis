import tkinter as tk
from tkinter import messagebox, filedialog
from DFWM.average_solutions import average_solutions
import re

class ParameterApp:
    def __init__(self, root):
        self.root = root
        root.title("Degenerated Four Wave Mixing")

        self.pretty_names = {
            "resolution": "rozdzielczość czasowa",
            "time_window": "okno czasowe [ps]",
            "z_saves": "rozdzielczość przestrzenna",
            "wavelength": "długość fali lasera [nm]",
            "fiber_length": "długość światłowodu [m]",
            "peak_power": "moc wiązki [W]",
            "rtol": "względna tolerancja (rtol)",
            "atol": "bezwzględna tolerancja (atol)",
            "n2": "współczynnik Kerra [m^2 / W]",
            "samples": "liczba modelowań",
            "neff_max": "maksymalny efektywny współczynnik załamania",
            "input_data_filepath": "ścieżka danych wejściowych"
        }

        self.params = {
            "resolution": 2**14,
            "time_window": 10,
            "z_saves": 51,
            "wavelength": 1064,
            "fiber_length": 0.15,
            "peak_power": 51000,
            "rtol": 1e-3,
            "atol": 1e-4,
            "n2": 2.6e-20,
            "samples": 1,
            "neff_max": 10,
            "input_data_filepath": "../data/neff_far_detuned_FWM.json"
        }

        self.entries = {}

        vcmd = (root.register(self.validate_numeric), "%P")

        row = 0
        for key, value in self.params.items():

            label_text = self.pretty_names.get(key, key)
            tk.Label(root, text=label_text).grid(
                row=row, column=0, sticky="w", padx=5, pady=3
            )

            entry = tk.Entry(root, width=40)

            if not isinstance(value, str):
                entry.config(validate="key", validatecommand=vcmd)

            entry.grid(row=row, column=1, padx=5, pady=3)
            entry.insert(0, str(value))
            self.entries[key] = entry

            if key == "input_data_filepath":
                browse_btn = tk.Button(root, text="Wybierz plik",
                                       command=self.browse_file)
                browse_btn.grid(row=row, column=2, padx=5, pady=3)
            row += 1

        self.start_button = tk.Button(root, text="Start", command=self.on_start)
        self.start_button.grid(row=row, column=0, columnspan=3, pady=10)

    def validate_numeric(self, new_value):
        if new_value == "":
            return True
        pattern = r"^-?(\d+(\.\d*)?|\.\d+)([eE]-?\d+)?$"
        return re.match(pattern, new_value) is not None

    def browse_file(self):
        filename = filedialog.askopenfilename(title="Wybierz plik")
        if filename:
            self.entries["input_data_filepath"].delete(0, tk.END)
            self.entries["input_data_filepath"].insert(0, filename)

    def on_start(self):
        self.toggle_widgets(state="disabled")

        params = {}
        for key, entry in self.entries.items():
            text = entry.get()
            default = self.params[key]

            if isinstance(default, str):
                params[key] = text
            else:
                try:
                    if "." in text or "e" in text.lower():
                        params[key] = float(text)
                    else:
                        params[key] = int(text)
                except ValueError:
                    params[key] = default

        try:
            average_solutions(**params)
        except Exception as e:
            messagebox.showerror("Error", f"An error occurred:\n{str(e)}")
        finally:
            self.toggle_widgets(state="normal")

    def toggle_widgets(self, state="disabled"):
        for e in self.entries.values():
            e.config(state=state)
        self.start_button.config(state=state)


if __name__ == "__main__":
    root = tk.Tk()
    app = ParameterApp(root)
    root.mainloop()
