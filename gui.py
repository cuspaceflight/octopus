"""
This tool can (eventually) be used to create, edit and save/load Plate objects,
and perform analysis of their mass flow distribution and (ambitiously)
combustion stability considerations from relationships found in the literature.

Perhaps it can be expanded to incorporate other Octopus features that would benefit from a GUI.

Need to check on changing coords from int to float
"""

from octopus import Fluid, Orifice, PropertySource, Manifold, Element
import tkinter as tk
from tkinter import ttk
from tkinter import messagebox

id = 0
id_array = []

fluids = {"Nitrous oxide": "nitrous oxide", "Isopropyl alcohol": "isopropyl alcohol"}
methods = {"Octopus Helmholtz (recommended)": "helmholz", "Python Thermo": "thermo"}
# As dicts in case we add fluids/methods that have awkward names

class Plate:
    """Class for representing the faceplate of an injector"""

    def __init__(self, id: int, diameter: float):
        """Initialise :class:`Plate` object.
        
        :param id: Plate ID (int)
        :param diameter: Plate diameter (m)
        
        """
        self.id = id
        self.diameter = diameter
        self.orifices = {"Z": []}
        self.orifice_count = 0
        self.x_centre = None
        self.y_centre = None

    def add_orifice(self, orifice: Orifice, x: int, y: int, arrayID="Z"):
        """Add a singular :class:'Orifice' object to the plate.
        Each orifice added has an array ID (string), if it was created as part of a series,
        and also a unique orifice ID (integer).
        The orifice's position is used for mass flow distribution calculations.

        :param Orifice: :class:'Orifice' to add to plate
        :param x: x position of Orifice on plate, relative to centre (int)
        :param y: y position of Orifice on plate, relative to centre (int)
        :param arrayID: Alphabetical identifier of orifice series (str, defaults to "Z")

        """
        self.orifice_count += 1

        if arrayID is None:
            self.orifices["Z"].append((orifice, self.orifice_count))
        else:
            if arrayID not in self.orifices:
                self.orifices[arrayID] = [(orifice, self.orifice_count)]
                print(f"New orifice array", arrayID)
                print(self.orifices[arrayID])
            else:
                self.orifices[arrayID].append((orifice, self.orifice_count))
                print(f"Orifice array {self.orifices[arrayID]} already exists")
                print(self.orifices[arrayID])

# Function for new injector plate, also has to update some elements of the GUI
# as this may be the first plate for this execution
def new_plate():
    try:
        diameter = float(diameter_entry.get())
        if diameter == 0:
            tk.messagebox.showerror("Error", "Plate diameter must be non-zero")
            return None
    except ValueError:
        tk.messagebox.showerror("Error", "Invalid type for plate diameter")
        return None

    global id, id_array, selected_plate
    id += 1
    id_array = list(range(1, id+1))

    plate_select_dropdown = tk.OptionMenu(frame_select_plate, selected_plate, *id_array)
    plate_select_dropdown.grid(row=0, column=1)
    plate_select_label["text"] = "Plate ID: "

    plate = Plate(id, diameter)

    plate_window = tk.Tk()
    plate_window.title(f"Plate ID {plate.id}")
    plate_window.columnconfigure(0, weight=1)
    plate_window.rowconfigure([0, 1], weight=1)

    plate_top_frame = tk.Frame(master=plate_window, height=30, padx=5, pady=5)
    plate_top_frame.grid(row=0, column=0)

    plate_parameters = tk.Label(master=plate_top_frame, height=3, text= \
        f"Injector ID: {id}\n Plate diameter: {plate.diameter} m")
    plate_parameters.pack()

    face_frame = tk.Frame(master=plate_window, width=300, height=300, padx=5, pady=5)
    face_frame.grid(row=1, column=0)

    face = tk.Canvas(master=face_frame, bg="white", height=600, width=600)
    face_outline = face.create_oval(50, 50, 550, 550, fill="black")
    Plate.x_centre, Plate.y_centre = 300, 300
    face.grid()
    
    plate_window.mainloop()

# Initialise the plate configuration window
window_cfg = tk.Tk()
window_cfg.iconphoto(True, tk.PhotoImage(file="img/favicon.png"))
window_cfg.title("Octopus injector pattern analysis")
window_cfg.resizable(False, False)
window_cfg.columnconfigure(0, weight=1)
window_cfg.rowconfigure([0, 1], weight=1)
window_cfg.rowconfigure([2, 3], weight=0)

# Tabs for different interfaces
tab_parent = ttk.Notebook(window_cfg)
tab_plate = ttk.Frame(tab_parent)
tab_manfiold = ttk.Frame(tab_parent)
tab_orifices = ttk.Frame(tab_parent)
tab_IO = ttk.Frame(tab_parent)

tab_parent.add(tab_plate, text="New plate")
tab_parent.add(tab_manfiold, text="New manifold")
tab_parent.add(tab_orifices, text="Add orifices")
tab_parent.add(tab_IO, text="Save/load")
tab_parent.grid(row=0, column=0, sticky="w")

# Details of the first tab, new plate
# Frame for title, possibly description in future
frame_plate_top = tk.Frame(master=tab_plate, width=50, height=5, padx=5, pady=5)
frame_plate_top.grid(row=1, column=0)

# Tab title
title_plate = tk.Label(master=frame_plate_top, text="Create a blank injector plate", \
                        borderwidth=2, relief="groove", height=1, width=50)
title_plate.config(font=(96))
title_plate.grid(row=0, column=0)

# Frame for plate diameter entry field and associated labels
frame_manage_diameter = tk.Frame(master=tab_plate, width=50, height=50, padx=5, pady=5)
frame_manage_diameter.grid(row=2, column=0)

diameter_entry = tk.Entry(master=frame_manage_diameter, width=10)
diameter_entry.grid(row=0, column=0, sticky="w")

diameter_entry_label = tk.Label(master=frame_manage_diameter, text="plate diameter, m")
diameter_entry_label.grid(row=0, column=1, sticky="e")

# Button to create a new plate with above parameters
new_plate_button = tk.Button(master=tab_plate, text="Create injector plate", command=new_plate, width=30, height=3)
new_plate_button.grid(row=3, column=0)

# Details of the third tab, add orifices
# Frame for title, possibly description in future
frame_orifice_top = tk.Frame(master=tab_orifices, width=50, height=5, padx=5, pady=5)
frame_orifice_top.grid(row=1, column=0)

# Tab title
title_orifices = tk.Label(master=frame_orifice_top, text="Add and remove orifices", \
                        borderwidth=2, relief="groove", height=1, width=50)
title_orifices.config(font=(96))
title_orifices.grid(row=0, column=0)

# Frame for plate selection option menu and associated labels
select_plate_frame = tk.Frame(master=tab_orifices, width=50, height=50, padx=5, pady=5)
select_plate_frame.grid(row=2, column=0)

# Option menu to choose from available plates
selected_plate = tk.IntVar()
plate_select_label = tk.Label(master=select_plate_frame, text="Plate ID: (Add or load a plate first)")
plate_select_label.grid(row=0, column=0)

# Frame for configuring a new orifice, which can be added one at a time or as part of a series
orifice_config_frame = tk.Frame(master=tab_orifices, width=50, height=5, padx=5, pady=5)
orifice_config_frame.grid(row=3, column = 0, sticky="w")

# Label for orifice type selection
orifice_config_label1 = tk.Label(master=orifice_config_frame, text="Orifice type:")
orifice_config_label1.grid(row=0, column=0, sticky="e")
orifice_config_label2 = tk.Label(master=orifice_config_frame, text="Orifice diameter:")
orifice_config_label2.grid(row=1, column=0, sticky="e")
orifice_config_label3 = tk.Label(master=orifice_config_frame, text="Manifold:")
orifice_config_label3.grid(row=2, column=0, sticky="e")

# Details of the second tab, manifold creation
# Frame for title, possibly description in future
frame_manifold_top = tk.Frame(master=tab_manfiold, width=50, height=5, padx=5, pady=5)
frame_manifold_top.grid(row=1, column=0)

# Tab title
title_manifold = tk.Label(master=frame_manifold_top, text="Create a new manifold", \
                          borderwidth=2, relief="groove", height=1, width=50)
title_manifold.config(font=(96))
title_manifold.grid(row=0, column=0)

# Frame for configuiring all manifold properties
manifold_frame = tk.Frame(master=tab_manfiold, width=50, height=5, padx=5, pady=5)
manifold_frame.grid(row=2, column = 0, sticky="w")

# Frame for selecting a fluid for the new manifold
fluid_select_frame = tk.Frame(master=manifold_frame, height=5, padx=5, pady=5)
fluid_select_frame.grid(row=0, column = 0, sticky="w")

# Fluid select label
manifold_label1 = tk.Label(master=fluid_select_frame, text="Manifold fluid:")
manifold_label1.grid(row=0, column=0, sticky="w")

# Fluid select menu
selected_fluid = tk.StringVar()
manifold_fluid_select = tk.OptionMenu(fluid_select_frame, selected_fluid, *fluids.keys())
manifold_fluid_select.grid(row=0, column=1, sticky="w")

# Frame for selecting a fluid property method
method_select_frame = tk.Frame(master=manifold_frame, height=5, padx=5, pady=5)
method_select_frame.grid(row=1, column=0, sticky="w")

# Fluid property method select label
manifold_label2 = tk.Label(master=method_select_frame, text="Fluid method:")
manifold_label2.grid(row=0, column=0, sticky="w")

# Fluid property method select menu
selected_method = tk.StringVar()
manifold_fluid_select = tk.OptionMenu(method_select_frame, selected_method, *methods.keys())
manifold_fluid_select.grid(row=0, column=1)

# Frame for setting fluid temperature in manifold
fluid_temp_frame = tk.Frame(master=manifold_frame, height=5, padx=5, pady=5)
fluid_temp_frame.grid(row=2, column=0, sticky="w")

# Manifold temperature label
manifold_label3 = tk.Label(master=fluid_temp_frame, text="Fluid temperature:")
manifold_label3.grid(row=0, column=0, sticky="w")

# Manifold temperature entry field
manifold_temp_entry = tk.Entry(master=fluid_temp_frame, width=10)
manifold_temp_entry.grid(row=0, column=1)

# Manifold temperature unit label
manifold_label4 = tk.Label(master=fluid_temp_frame, text="K")
manifold_label4.grid(row=0, column=3)

# Frame for setting fluid pressure in manifold
fluid_pressure_frame = tk.Frame(master=manifold_frame, height=5, padx=5, pady=5)
fluid_pressure_frame.grid(row=3, column=0, sticky="w")

# Manifold pressure label
manifold_label5 = tk.Label(master=fluid_pressure_frame, text="Fluid pressure:")
manifold_label5.grid(row=0, column=0, sticky="w")

# Manifold pressure entry field
manifold_pressure_entry = tk.Entry(master=fluid_pressure_frame, width=10)
manifold_pressure_entry.grid(row=0, column=1)

# Manifold pressure unit label
manifold_label6 = tk.Label(master=fluid_pressure_frame, text="bar")
manifold_label6.grid(row=0, column=3)

window_cfg.mainloop()