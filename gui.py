"""
This tool can (eventually) be used to create, edit and save/load Plate objects,
and perform analysis of their mass flow distribution and (ambitiously)
combustion stability considerations from relationships found in the literature.

Perhaps it can be expanded to incorporate other Octopus features that would benefit from a GUI.

Need to check on changing coords from int to float
"""

from octopus import Fluid, Orifice, PropertySource, Manifold, Element
import tkinter as tk
import numpy as np
from tkinter import ttk
from tkinter import messagebox

# Global variables
last_plate_id = 0
Omodel = None
Otype = None
Oseries = None
Mfluid = None
Mmethod = None
plate_id_list = []
plates = [0] # Initial 0 is so plates[plate_id] behaves intuitively
last_manifold_id = 0
manifold_id_list = []
manifolds = [0] # Initial 0 is so plates[plate_id] behaves intuitively
manifold_mdot_list = [0] # Initial 0 is so manifold_mdot_list[manifold_id] behaves intuitively

class Plate:
    """Class for representing the faceplate of an injector"""

    def __init__(self, id: int, diameter: float, pc: float, diameter_px: int):
        """Initialise :class:`Plate` object.
        
        :param id: Plate ID (int)
        :param diameter: Plate diameter (m)
        :param pc: Combustion chamber pressure (Pa)
        :param diameter_x: Injector plate display diameter (pixels)
        
        """
        self.id = id
        self.diameter = diameter
        self.pc = pc
        self.orifices = {"Z": []}
        self.orifice_count = 0
        self.x_centre = None
        self.y_centre = None
        self.diameter_px = diameter_px

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

    def add_window(self):
        # Create the display window for this plate
        self.window = tk.Tk()
        self.window.resizable(False, False)
        self.window.title(f"Plate ID {self.id}")

        # Frame for the plate title and parameter list
        self.plate_top_frame = tk.Frame(master=self.window, height=30, padx=5, pady=5)
        self.plate_top_frame.grid(row=0, column=0)

        # Label for the parameters
        self.plate_parameters = tk.Label(master=self.plate_top_frame, height=3, text= \
            f"Plate ID: {self.id}\n Plate diameter: {self.diameter} m\n"
            f"Chamber pressure: {self.pc/1E5} bar")
        self.plate_parameters.grid(row=0, column=0)

        # Frame for the plate canvas
        self.face_frame = tk.Frame(master=self.window, width=300, height=300, padx=5, pady=5)
        self.face_frame.grid(row=1, column=0)

        # Create the canvas for representing the plate and its orifices
        self.face = tk.Canvas(master=self.face_frame, bg="white", \
                              height=self.diameter_px+10, width=self.diameter_px+10)
        self.face_outline = self.face.create_oval(6, 6, self.diameter_px+6, self.diameter_px+6, fill="black")
        self.x_centre, self.y_centre = (self.diameter_px+10)/2, (self.diameter_px+10)/2

        # Sneakily create a preview at the centre so it can be updated later, without creating a new object
        # hence there is only ever one preview and we don't need to track if one has been created
        # (not sure how best to do this for a series)
        self.preview = self.face.create_oval(self.x_centre, self.y_centre, \
                                             self.x_centre, self.y_centre, fill="white")
        self.face.itemconfigure(self.preview, state="hidden")
        self.face.grid()
        
        self.window.mainloop()


# Function for new injector plate, also has to update some elements of the GUI
# as this may be the first plate for this execution
def new_plate():
    # Attempt to retrieve the needed parameters
    try:
        diameter = float(diameter_entry.get())
        pc = float(chamber_pressure_entry.get())
        diameter_px = float(canvas_size_entry.get())
        if diameter <= 0 or pc <= 0:
            tk.messagebox.showerror("New plate error", "Plate diameter and chamber pressure must be greater than zero")
            return None
    except Exception:
        tk.messagebox.showerror("New plate error", "Invalid plate parameters")
        return None

    global last_plate_id, plate_id_list, selected_plate, plate_select_dropdown

    # Increment the ID value and update the list of IDs
    last_plate_id += 1
    plate_id_list = list(range(1, last_plate_id+1))

    # Update the option menu in the add orifice tab with the new ID list and
    # remove the warning that there are no plates - one has just been added
    plate_select_dropdown = tk.OptionMenu(select_plate_frame, selected_plate, *plate_id_list)
    plate_select_dropdown.grid(row=0, column=1)
    plate_select_label["text"] = "Plate ID: "

    # Create the plate object, append it to the list
    # for access by indexing with plate id
    # (should change this to use a dict)
    plates.append(Plate(last_plate_id, diameter, pc*1E5, diameter_px))

    # Create and display the plate's window
    # This is handled within the object, so other methods can reference
    # specific plate IDs appropriately
    plates[last_plate_id].add_window()


# Function for new manifold, like the above has to update
# some elements of the GUI
def new_manifold():
    global Mfluid, Mmethod
    try:
        fluid_id = Mfluid
        method_id = Mmethod
        p = float(manifold_pressure_entry.get())*1E5
        T = float(manifold_temp_entry.get())
        if p == 0 or T == 0:
            tk.messagebox.showerror("New manifold error", "Manifold temperature and pressure must be greater than 0")
            return None
    except Exception:
        tk.messagebox.showerror("New manifold error", "Invalid manifold parameters")
        return None

    global last_manifold_id, manifold_id_list, selected_manifold, manifold_mdot_list, manifold_select_dropdown

    # Increment the ID value and update the list of IDs
    last_manifold_id += 1
    id = last_manifold_id
    manifold_id_list = list(range(1, last_manifold_id+1))

    # Update the option menu in the add orifice tab with the new ID list and
    # remove the warning that there are no plates - one has just been added
    manifold_select_dropdown = tk.OptionMenu(select_manifold_frame, selected_manifold, *manifold_id_list)
    manifold_select_dropdown.grid(row=0, column=1)
    manifold_select_label["text"] = "Manifold ID: "

    # Create the fluid and manifold objects
    fluid = Fluid(ID=fluid_id, T=T, P=p, method=method_id)
    manifold = Manifold(fluid, PropertySource(p=p, T=T))
    manifolds.append(manifold)
    manifold_mdot_list.append(0)

    # Create the window for this manifold
    manifold_window = tk.Tk()
    manifold_window.resizable(False, False)
    manifold_window.title(f"Manifold ID {id}")

    # Frame for the plate title and parameter list
    manifold_top_frame = tk.Frame(master=manifold_window, height=30, padx=5, pady=5)
    manifold_top_frame.grid(row=0, column=0)

    # Label for the fixed parameters
    plate_static_label = tk.Label(master=manifold_top_frame, height=4, text= \
        f"Manifold ID: {id}\n Manifold fluid: {manifold.fluid.ID}\n"
        f"Manifold pressure: {manifold.p/1E5} bar\n Manifold temperature: {manifold.T} K")
    plate_static_label.grid(row=0, column=0, sticky="w")
    
    # Label for the masss flow rate, updated when orifices are added
    plate_massflow_label = tk.Label(master=manifold_top_frame, text= \
        f"Manifold mass flow rate: {manifold_mdot_list[id]} kg/s")
    plate_massflow_label.grid(row=1, column=0, sticky="w")

    manifold_window.mainloop()


# Update the mass flow model when a radio button is clicked in the orifice tab
def model_update():
    global Omodel
    Omodel = selected_model.get()


# Update the orifice type when a radio button is clicked in the orifice tab
# If the type is not Waxman, disable that mass flow model
# If the type is Waxman, select that as the mass flow model and disable
# the other mass flow models. At some point, should add the other models for
# a Waxman injector, using the throat as the diameter
def type_update():
    global Otype, Omodel
    Otype = orifice_type.get()

    if Otype == 0: # (Straight orifice)
        orifice_model_WAXMAN.config(state="disabled")
        orifice_model_SPI.config(state="normal")
        orifice_model_SPI.select()
        Omodel = "SPI"
        orifice_model_HEM.config(state="normal")
        orifice_model_DYER.config(state="normal")
    
    if Otype == 1: # (cavitating orifice)
        orifice_model_WAXMAN.config(state="normal")
        orifice_model_WAXMAN.select()
        Omodel = "WAXMAN"
        orifice_model_SPI.config(state="disabled")
        orifice_model_HEM.config(state="disabled")
        orifice_model_DYER.config(state="disabled")


# Update the fluid when a radio button is clicked in the manifold tab
# If the fluid is IPA, disable the Octopus Helmholtz EOS option
def fluid_update():
    global Mfluid, Mmethod
    Mfluid = selected_fluid.get()

    if Mfluid == "nitrous oxide":
        manifold_method_Helmholtz.select()
        Mmethod = "helmholz"
        manifold_method_Helmholtz.config(state="normal")

    if Mfluid == "isopropyl alcohol":
        manifold_method_thermo.select()
        Mmethod = "thermo"
        manifold_method_Helmholtz.config(state="disabled")


def method_update():
    global Mmethod
    Mmethod = selected_method.get()


# Disable all the orifice configuration options once the settings are
# confirmed. This makes it easier to draw the plate preview.
# This also enables all the orifice creation options, updates preview
def orifice_confirm():
    global plate_select_dropdown, manifold_select_dropdown
    # Slightly clumsy error handling

    try:
        float(orifice_diameter_entry.get())
        if float(orifice_diameter_entry.get()) <= 0:
            tk.messagebox.showerror("Orifice settings error", "Orifice diameter must be greater than zero")
            return None
    except Exception:
        tk.messagebox.showerror("Orifice settings error", "Invalid orifice diameter")
        return None

    try:
        int(selected_plate.get())
        if selected_plate.get() == 0:
            tk.messagebox.showerror("Orifice settings error", "No plate selected")
            return None
    except Exception:
        tk.messagebox.showerror("Orifice settings error", "No plate selected")
        return None

    try:
        int(selected_manifold.get())
        if selected_manifold.get() == 0:
            tk.messagebox.showerror("Orifice settings error", "No manifold selected")
            return None
    except Exception:
        tk.messagebox.showerror("Orifice settings error", "No manifold selected")
        return None

    orifice_type_straight.config(state="disabled")
    orifice_type_waxman.config(state="disabled")

    orifice_diameter_entry.config(state="disabled")

    orifice_model_SPI.config(state="disabled")
    orifice_model_HEM.config(state="disabled")
    orifice_model_DYER.config(state="disabled")
    orifice_model_WAXMAN.config(state="disabled")

    plate_select_dropdown.config(state="disabled")
    manifold_select_dropdown.config(state="disabled")

    """orifice_add_single.config(state="normal")
    orifice_add_series.config(state="normal")"""
    orifice_count_entry.config(state="normal")

    orifice_config_edit.config(state="normal")
    orifice_config_confirm.config(state="disabled")

    new_orifice_aentry.config(state="normal")
    new_orifice_rentry.config(state="normal")

    tab_parent.tab(0, state="disabled")
    tab_parent.tab(1, state="disabled")
    tab_parent.tab(3, state="disabled")

    orifice_preview_update()


# Enable all the orifice configuration options again if the user
# wants to abort creating a new orifice. This disables all the
# orifice creation options again, and clears the preview
def orifice_edit():
    global plate_select_dropdown, manifold_select_dropdown

    type_update()
    orifice_type_straight.config(state="normal")
    orifice_type_waxman.config(state="normal")

    orifice_diameter_entry.config(state="normal")

    model_update()

    plate_select_dropdown.config(state="normal")
    manifold_select_dropdown.config(state="normal")

    """orifice_add_single.config(state="disabled")
    orifice_add_series.config(state="disabled")"""
    orifice_count_entry.config(state="disabled")

    orifice_config_confirm.config(state="normal")
    orifice_config_edit.config(state="disabled")

    new_orifice_aentry.config(state="disabled")
    new_orifice_rentry.config(state="disabled")

    # Find the old plate ID, hide its preview
    plate = plates[selected_plate.get()]
    plate.face.itemconfigure(plate.preview, state="hidden")

    tab_parent.tab(0, state="normal")
    tab_parent.tab(1, state="normal")
    tab_parent.tab(3, state="normal")


# Update menus for either a single orifice or series of orifices
# in the orifice tab when a radio button is clicked
def single_series_update():
    global Oseries
    Oseries = orifice_single_series.get()

    if Oseries == "single":
        pass

    if Oseries == "series":
        pass


# Whenever the x position field is updated, try and update
# the canvas preview.
def orifice_preview_update(*args):
    try:
        angle = float(new_orifice_ang.get())
        radius = float(new_orifice_r.get())
        od = float(orifice_diameter_entry.get())
        x = np.sin(angle)*radius
        y = np.cos(angle)*radius
    except Exception:
        return None

    # Polar to Cartesian for the first orifice
    x = np.sin(angle*np.pi/180)*radius
    y = np.cos(angle*np.pi/180)*radius

    # Get the plate we're working with
    plate_id = int(selected_plate.get())
    plate = plates[plate_id]

    # And the manifold
    manifold_id = int(selected_manifold.get())
    manifold = manifolds[manifold_id]
    
    # For IPA (fuel), set red preview colour, for N2O (ox), set blue
    # else leave it white
    if manifold.fluid.ID == "nitrous oxide":
        plate.face.itemconfig(plate.preview, fill="blue", state="normal")
    elif manifold.fluid.ID == "isopropyl alcohol":
        plate.face.itemconfig(plate.preview, fill="red", state="normal")
    else:
        plate.face.itemconfig(plate.preview, fill="white", state="normal")


    # Show the preview
    plate.face.itemconfigure(plate.preview, state="normal")

    # Find the mm per pixel
    mmpp = 1000*plate.diameter/plate.diameter_px
    # Find the orifice diameter in pixels
    odpx = od/mmpp
    
    # Update the preview, with appropriately scaled size
    plate.face.coords(plate.preview, plate.x_centre+x, plate.y_centre-y, plate.x_centre+x-odpx, plate.y_centre-y-odpx)
    plate.face.grid()


# Initialise the configuration window
window_cfg = tk.Tk()
window_cfg.iconphoto(True, tk.PhotoImage(file="img/favicon.png"))
window_cfg.title("Octopus injector design")
window_cfg.resizable(False, False)

# Tabs for different interfaces
tab_parent = ttk.Notebook(window_cfg)
tab_plate = ttk.Frame(tab_parent)
tab_manifold = ttk.Frame(tab_parent)
tab_orifices = ttk.Frame(tab_parent)
tab_IO = ttk.Frame(tab_parent)

tab_parent.add(tab_plate, text="New plate")
tab_parent.add(tab_manifold, text="New manifold")
tab_parent.add(tab_orifices, text="Manage orifices")
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
manage_diameter_frame = tk.Frame(master=tab_plate, width=50, height=50, padx=5, pady=5)
manage_diameter_frame.grid(row=2, column=0, sticky="w")

# Plate diameter entry field label
diameter_entry_label = tk.Label(master=manage_diameter_frame, text="Plate diameter:")
diameter_entry_label.grid(row=0, column=0, sticky="w")

# Plate diameter entry field
diameter_entry = tk.Entry(master=manage_diameter_frame, width=5)
diameter_entry.grid(row=0, column=1, sticky="w")

# Label for plate diameter unit
diameter_entry_unit = tk.Label(master=manage_diameter_frame, text="m")
diameter_entry_unit.grid(row=0, column=2, sticky="w")

# Frame for combustion chamber pressure entry field and associated labels
chamber_pressure_frame = tk.Frame(master=tab_plate, width=50, height=50, padx=5, pady=5)
chamber_pressure_frame.grid(row=3, column=0, sticky="w")

# Chamber pressure entry field label
chamber_pressure_label = tk.Label(master=chamber_pressure_frame, text="Combustion chamber pressure:")
chamber_pressure_label.grid(row=0, column=0, sticky="w")

# Chamber pressure entry field
chamber_pressure_entry = tk.Entry(master=chamber_pressure_frame, width=5)
chamber_pressure_entry.grid(row=0, column=1, sticky="w")

# Label for chamber pressure unit
chamber_pressure_unit = tk.Label(master=chamber_pressure_frame, text="bar")
chamber_pressure_unit.grid(row=0, column=2, sticky="w")

# Frame for setting canvas size
canvas_size_frame = tk.Frame(master=tab_plate, width=50, height=50, padx=5, pady=5)
canvas_size_frame.grid(row=4, column=0, sticky="w")

# Canvas size label
canvas_size_label = tk.Label(master=canvas_size_frame, text="Injector display diameter:")
canvas_size_label.grid(row=0, column=0, sticky="w")

# Canvas size entry field
canvas_size_entry = tk.Entry(master=canvas_size_frame, width=5)
canvas_size_entry.grid(row=0, column=1, sticky="w")

# Canvas size unit label
canvas_size_unit = tk.Label(master=canvas_size_frame, text="px")
canvas_size_unit.grid(row=0, column=2, sticky="w")

# Button to create a new plate with above parameters (opens in new window)
new_plate_button = tk.Button(master=tab_plate, text="Create injector plate", command=new_plate, width=30, height=3)
new_plate_button.grid(row=5, column=0)


# Details of the third tab, add orifices
# Frame for title, possibly description in future
frame_orifice_top = tk.Frame(master=tab_orifices, width=50, height=5, padx=5, pady=5)
frame_orifice_top.grid(row=1, column=0)

# Tab title
title_orifices = tk.Label(master=frame_orifice_top, text="Add and remove orifices", \
                        borderwidth=2, relief="groove", height=1, width=50)
title_orifices.config(font=(96))
title_orifices.grid(row=0, column=0)

# Frame for configuring a new orifice
orifice_config_frame = tk.Frame(master=tab_orifices, width=50, height=5, padx=5, pady=5,\
                                relief="sunken", bd=2)
orifice_config_frame.grid(row=3, column = 0, sticky="w")

# Frame for selection of orifice type
orifice_type_frame = tk.Frame(master=orifice_config_frame, width=50, height=50, padx=5, pady=5)
orifice_type_frame.grid(row=0, column=0, sticky="w")

# Label for orifice type selection
orifice_type_label = tk.Label(master=orifice_type_frame, text="Orifice type:")
orifice_type_label.grid(row=0, column=0, sticky="w")

# Radio buttons for orifice type selection
orifice_type = tk.IntVar(value=0)
orifice_type_straight = tk.Radiobutton(master=orifice_type_frame, variable=orifice_type,\
                                       text="Classic (straight)", value=0, command=type_update)
orifice_type_straight.grid(row=1, column=0, sticky="w")

orifice_type_waxman = tk.Radiobutton(master=orifice_type_frame, variable=orifice_type,\
                                       text="Cavitating (Waxman) ", value=1, command=type_update)
orifice_type_waxman.grid(row=1, column=1, sticky="w")

# Frame for orifice diameter entry field and associated labels
orifice_diameter_frame = tk.Frame(master=orifice_config_frame, width=50, height=50, padx=5, pady=5)
orifice_diameter_frame.grid(row=1, column=0, sticky="w")

# Label for orifice diameter
orifice_diameter_label = tk.Label(master=orifice_diameter_frame, text="Orifice exit diameter:")
orifice_diameter_label.grid(row=0, column=0, sticky="w")

# Entry field for orifice diameter
orifice_diameter_entry = tk.Entry(master=orifice_diameter_frame, width=10)
orifice_diameter_entry.grid(row=0, column=1, sticky="w")

# Label for orifice diameter unit
orifice_diameter_unit = tk.Label(master=orifice_diameter_frame, text="mm")
orifice_diameter_unit.grid(row=0, column=2, sticky="w")

# Frame for mass flow rate model selection
orifice_model_frame = tk.Frame(master=orifice_config_frame, width=50, height=50, padx=5, pady=5)
orifice_model_frame.grid(row=2, column=0, sticky="w")

# Label for mass flow rate model selection
orifice_model_label = tk.Label(master=orifice_model_frame, text="Mass flow model:")
orifice_model_label.grid(row=0, column=0, sticky="w")

# Radio buttons for mass flow model selection
selected_model = tk.StringVar(value="SPI")
orifice_model_SPI = tk.Radiobutton(master=orifice_model_frame, variable=selected_model,\
                                   text="Single phase incompressible",\
                                   value="SPI", command=model_update)
orifice_model_SPI.grid(row=1, column=0, sticky="w")

orifice_model_HEM = tk.Radiobutton(master=orifice_model_frame, variable=selected_model,\
                                   text="Homogenous equlibrium",\
                                   value="HEM", command=model_update)
orifice_model_HEM.grid(row=1, column=1, sticky="w")

orifice_model_DYER = tk.Radiobutton(master=orifice_model_frame, variable=selected_model,\
                                   text="Dyer (Solomon corrected)",\
                                   value="DYER", command=model_update)
orifice_model_DYER.grid(row=2, column=0, sticky="w")

orifice_model_WAXMAN = tk.Radiobutton(master=orifice_model_frame, variable=selected_model,\
                                      text="Cavitating (Waxman)",\
                                      value="WAXMAN", command=model_update)
orifice_model_WAXMAN.grid(row=2, column=1, sticky="w")
# Call model_update and type_update to setup with default values
model_update()
type_update()

# Frame for plate selection option menu and associated labels
select_plate_frame = tk.Frame(master=orifice_config_frame, width=50, height=50, padx=5, pady=5)
select_plate_frame.grid(row=3, column=0, sticky="w")

# Option menu to choose from available plates. OptionMenu widget is added / updated
# by new_plate() / import_manifold()
selected_plate = tk.IntVar()
plate_select_label = tk.Label(master=select_plate_frame, text="Plate ID: (Add or load a plate first)")
plate_select_label.grid(row=0, column=0, sticky="w")

# Frame for selection of orifice parent manifold
select_manifold_frame = tk.Frame(master=orifice_config_frame, width=50, height=50, padx=5, pady=5)
select_manifold_frame.grid(row=4, column=0, sticky="w")

# Option menu to choose from available manifolds. OptionMenu widget is added / updated
# by new_manifold() / import_manifold()
selected_manifold = tk.StringVar()
manifold_select_label = tk.Label(master=select_manifold_frame, text="Manifold ID: (Add or load a manifold first)")
manifold_select_label.grid(row=0, column=0, sticky="w")

# Frame for buttons to confirm / edit orifice settings
orifice_config_buttons = tk.Frame(master=orifice_config_frame, width=50, height=50, padx=5, pady=5)
orifice_config_buttons.grid(row=5, column=0)

# Buttons for confirming, editing orifice settings
orifice_config_confirm = tk.Button(master=orifice_config_buttons, text="Confirm settings", command=orifice_confirm)
orifice_config_confirm.grid(row=0, column=0)
orifice_config_edit = tk.Button(master=orifice_config_buttons, text="Edit settings", command=orifice_edit)
orifice_config_edit.grid(row=0, column=1)
# Disable this initially as config need to be confirmed first
orifice_config_edit.config(state="disabled")

# Frame for orifice creation options
orifice_create_frame = tk.Frame(master=tab_orifices, padx=5, pady=5, relief="sunken", bd=2)
orifice_create_frame.grid(row=6, column=0, sticky="w")

"""
# Frame for choosing whether to add an orifice as a series or on its own
orifice_add_type = tk.Frame(master=orifice_create_frame, padx=5, pady=5)
orifice_add_type.grid(row=0, column=0, sticky="w")

# Label for single / series orifice
orifice_single_label = tk.Label(master=orifice_add_type, text="Orifices to add:")
orifice_single_label.grid(row=0, column=0, sticky="w")

# Radio buttons for new orifice individual / series selection
orifice_single_series = tk.StringVar(value="single")
orifice_add_single = tk.Radiobutton(master=orifice_add_type, variable=orifice_single_series,\
                                   text="Single orifice",\
                                   value="single", command=single_series_update)
orifice_add_single.grid(row=1, column=0, sticky="w")

orifice_add_series = tk.Radiobutton(master=orifice_add_type, variable=orifice_single_series,\
                                   text="Series of orifices",\
                                   value="pattern", command=single_series_update)
orifice_add_series.grid(row=1, column=1, sticky="w")
# Call these to setup default states / values
single_series_update()
orifice_add_single.config(state="disabled")
orifice_add_series.config(state="disabled")
"""

# Frame for choosing number of orifices to add
orifice_count_frame = tk.Frame(master=orifice_create_frame, padx=5, pady=5)
orifice_count_frame.grid(row=0, column=0, sticky="w")

# Label for orifice count
orifice_count_label = tk.Label(master=orifice_count_frame, text="Number of new orifices:")
orifice_count_label.grid(row=0, column=0, sticky="w")

# Entry field for orifice count
orifice_new_count = tk.StringVar()
orifice_count_entry = tk.Entry(master=orifice_count_frame, textvariable=orifice_new_count, width=5)
orifice_count_entry.grid(row=0, column=1, sticky="w")
# Disable initially
orifice_count_entry.config(state="disabled")


# Frame for orifice position (for first in series if a series is elected)
orifice_position_frame = tk.Frame(master=orifice_create_frame, padx=5, pady=5)
orifice_position_frame.grid(row=7, column=0, sticky="w")

# Label for orifice position description
orifice_position_label = tk.Label(master=orifice_position_frame, text="Polar position of first orifice:")
orifice_position_label.grid(row=0, column=0, sticky="w")

# angle entry field label
orifice_position_alabel = tk.Label(master=orifice_position_frame, text="Angle:")
orifice_position_alabel.grid(row=0, column=1, sticky="w")

# angle entry field
# Use a DoubleVar so we know when this is updated
new_orifice_ang = tk.DoubleVar()
new_orifice_ang.trace_add("write", orifice_preview_update)
new_orifice_aentry = tk.Entry(master=orifice_position_frame, textvariable=new_orifice_ang, width=5)
new_orifice_aentry.grid(row=0, column=2, sticky="w")

# Angle entry unit label
orifice_aentry_unit = tk.Label(master=orifice_position_frame, text=u"\N{DEGREE SIGN}")
orifice_aentry_unit.grid(row=0, column=3, sticky="w")

# Radius entry field label
orifice_position_rrlabel = tk.Label(master=orifice_position_frame, text="Radius:")
orifice_position_rrlabel.grid(row=0, column=4, sticky="w")

# Radius entry field
# Use a DoubleVar so we know when this is updated
new_orifice_r = tk.DoubleVar()
new_orifice_r.trace("w", orifice_preview_update)
new_orifice_rentry = tk.Entry(master=orifice_position_frame, textvariable=new_orifice_r, width=5)
new_orifice_rentry.grid(row=0, column=5, sticky="w")

# Radius entry unit label
orifice_rentry_unit = tk.Label(master=orifice_position_frame, text="mm")
orifice_rentry_unit.grid(row=0, column=6, sticky="w")

# Set up default states
new_orifice_aentry.config(state="disabled")
new_orifice_rentry.config(state="disabled")

# Details of the second tab, manifold creation
# Frame for title, possibly description in future
frame_manifold_top = tk.Frame(master=tab_manifold, width=50, height=5, padx=5, pady=5)
frame_manifold_top.grid(row=1, column=0)

# Tab title
title_manifold = tk.Label(master=frame_manifold_top, text="Create a new manifold", \
                        borderwidth=2, relief="groove", height=1, width=50)
title_manifold.config(font=(96))
title_manifold.grid(row=0, column=0)

# Frame for configuiring all manifold properties
manifold_frame = tk.Frame(master=tab_manifold, width=50, height=5, padx=5, pady=5)
manifold_frame.grid(row=2, column = 0, sticky="w")

# Frame for selecting a fluid for the new manifold
fluid_select_frame = tk.Frame(master=manifold_frame, height=5, padx=5, pady=5)
fluid_select_frame.grid(row=0, column = 0, sticky="w")

# Fluid select label
manifold_fluid_label = tk.Label(master=fluid_select_frame, text="Manifold fluid:")
manifold_fluid_label.grid(row=0, column=0, sticky="w")

# Fluid select radio buttons
selected_fluid = tk.StringVar(value="nitrous oxide")
manifold_fluid_N2O = tk.Radiobutton(master=fluid_select_frame, variable=selected_fluid,\
                                    text="Nitrous oxide", value="nitrous oxide",\
                                    command=fluid_update)
manifold_fluid_N2O.grid(row=1, column=0, sticky="w")
manifold_fluid_IPA = tk.Radiobutton(master=fluid_select_frame, variable=selected_fluid,\
                                    text="Isopropyl alcohol", value="isopropyl alcohol",\
                                    command=fluid_update)
manifold_fluid_IPA.grid(row=1, column=1, sticky="w")

# Frame for selecting a fluid property method
method_select_frame = tk.Frame(master=manifold_frame, height=5, padx=5, pady=5)
method_select_frame.grid(row=1, column=0, sticky="w")

# Fluid property method select label
manifold_method_label = tk.Label(master=method_select_frame, text="Fluid method:")
manifold_method_label.grid(row=0, column=0, sticky="w")

# Fluid property method radio buttons
selected_method = tk.StringVar(value="helmholz")
manifold_method_thermo = tk.Radiobutton(master=method_select_frame, variable=selected_method,\
                                        text="Python Thermo", value="thermo",\
                                        command=method_update)
manifold_method_thermo.grid(row=1, column=0, sticky="w")
manifold_method_Helmholtz = tk.Radiobutton(master=method_select_frame, variable=selected_method,\
                                          text="Octopus Helmholtz", value="helmholz",\
                                          command=method_update)
manifold_method_Helmholtz.grid(row=1, column=1, sticky="w")
# Call method_update and fluid_update to setup default values
method_update()
fluid_update()

# Frame for setting fluid temperature in manifold
fluid_temp_frame = tk.Frame(master=manifold_frame, height=5, padx=5, pady=5)
fluid_temp_frame.grid(row=2, column=0, sticky="w")

# Manifold temperature label
manifold_temp_label = tk.Label(master=fluid_temp_frame, text="Fluid temperature:")
manifold_temp_label.grid(row=0, column=0, sticky="w")

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
manifold_pressure_label = tk.Label(master=fluid_pressure_frame, text="Fluid pressure:")
manifold_pressure_label.grid(row=0, column=0, sticky="w")

# Manifold pressure entry field
manifold_pressure_entry = tk.Entry(master=fluid_pressure_frame, width=10)
manifold_pressure_entry.grid(row=0, column=1)

# Manifold pressure unit label
manifold_pressure_unit = tk.Label(master=fluid_pressure_frame, text="bar")
manifold_pressure_unit.grid(row=0, column=3)

# Button to create a new manifold with above parameters (opens in new window)
manifold_button = tk.Button(master=tab_manifold, text="Create manifold", command=new_manifold, width=30, height=3)
manifold_button.grid(row=3, column=0)


window_cfg.mainloop()
