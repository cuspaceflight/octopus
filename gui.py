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
                print(f"New orifice array", arrayID)
                self.orifices[arrayID] = [(orifice, self.orifice_count)]
                print(self.orifices[arrayID])
            else:
                self.orifices[arrayID].append((orifice, self.orifice_count))
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

    plate_select_dropdown = tk.OptionMenu(frame_select_plate, selected_plate, 1, *id_array[1:])
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

# Tabs for creating a new plate, adding orificies, saving and loading plates
tab_parent = ttk.Notebook(window_cfg)
tab_manage = ttk.Frame(tab_parent)
tab_orifices = ttk.Frame(tab_parent)
tab_IO = ttk.Frame(tab_parent)

tab_parent.add(tab_manage, text="New plate")
tab_parent.add(tab_orifices, text="Add orifices")
tab_parent.add(tab_IO, text="Import/export")
tab_parent.grid(row=0, column=0, sticky="w")

# Details of the first tab, new plate
frame_manage_top = tk.Frame(master=tab_manage, width=50, height=5, padx=5, pady=5)
frame_manage_top.grid(row=1, column=0)

title_manage = tk.Label(master=frame_manage_top, text="Create a blank injector plate", \
                        borderwidth=2, relief="groove", height=1, width=50)
title_manage.config(font=(96))
title_manage.grid(row=0, column=0)

frame_manage_diameter = tk.Frame(master=tab_manage, width=50, height=50, padx=5, pady=5)
frame_manage_diameter.grid(row=2, column=0)

diameter_entry = tk.Entry(master=frame_manage_diameter, width=10)
diameter_entry.grid(row=0, column=0, sticky="w")

diameter_entry_label = tk.Label(master=frame_manage_diameter, text="plate diameter, m")
diameter_entry_label.grid(row=0, column=1, sticky="e")

new_plate_button = tk.Button(master=tab_manage, text="Create injector plate", command=new_plate, width=30, height=3)
new_plate_button.grid(row=3, column=0)

# Details of the second tab, add orifices
# Example config pulled from /examples/dev.py
nitrous = Fluid('nitrous oxide', P=20e5, method='helmholz')
isopropanol = Fluid('isopropanol', P=18e5, T=400, method='thermo')
nitrous_manifold = Manifold(nitrous, PropertySource(p=18e5, T=250))
ipa_manifold = Manifold(isopropanol, PropertySource(p=18e5, T=400))
nitrous_orifice = Orifice(nitrous_manifold, 1e-2, 2e-3)
ipa_orifice = Orifice(ipa_manifold, 1e-2, 1e-3)
element = Element([nitrous_orifice, nitrous_orifice], [ipa_orifice, ipa_orifice])

frame_orifice_top = tk.Frame(master=tab_orifices, width=50, height=5, padx=5, pady=5)
frame_orifice_top.grid(row=1, column=0)

title_orifices = tk.Label(master=frame_orifice_top, text="Add and remove orifices", \
                        borderwidth=2, relief="groove", height=1, width=50)
title_orifices.config(font=(96))
title_orifices.grid(row=0, column=0)

frame_select_plate = tk.Frame(master=tab_orifices, width=50, height=50, padx=5, pady=5)
frame_select_plate.grid(row=2, column=0)

selected_plate = tk.IntVar(frame_select_plate)
plate_select_label = tk.Label(master=frame_select_plate, text="Plate ID: (Add or load a plate first)")
plate_select_label.grid(row=0, column=0)

single_orifice_frame = tk.Frame(master=tab_orifices, width=50, height=5, padx=5, pady=5)
single_orifice_frame.grid(row=3, column = 0, sticky="w")

single_orifice_label1 = tk.Label(master=single_orifice_frame, text="Orifice type:")
single_orifice_label1.grid(row=0, column=0, sticky="e")
single_orifice_label2 = tk.Label(master=single_orifice_frame, text="Orifice diameter:")
single_orifice_label2.grid(row=1, column=0, sticky="e")
single_orifice_label3 = tk.Label(master=single_orifice_frame, text="Manifold:")
single_orifice_label3.grid(row=2, column=0, sticky="e")


window_cfg.mainloop()