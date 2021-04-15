from octopus import Fluid, Orifice, PropertySource, Manifold, Element
import tkinter as tk
from tkinter import ttk

id = 0

def new_plate():
    try:
        diameter = float(diameter_entry.get())
        if diameter == 0:
            print("Invalid diameter for injector plate")
            return None
    except ValueError:
        print("Invalid diameter for injector plate")
        return None

    global id
    id += 1

    plate_window = tk.Tk()
    plate_window.title(f"Plate ID {id}")
    plate_window.columnconfigure(0, weight=1)
    plate_window.rowconfigure([0, 1], weight=1)

    plate_top_frame = tk.Frame(master=plate_window, height=30, padx=5, pady=5)
    plate_top_frame.grid(row=0, column=0)

    plate_title = tk.Label(master=plate_top_frame, text="Injector plate configuration")
    plate_title.pack()

    plate_parameters = tk.Label(master=plate_top_frame, height=3, text= \
        f"Injector ID: {id}\n Plate diameter: {diameter} mm")
    plate_parameters.pack()

    face_frame = tk.Frame(master=plate_window, width=300, height=300, padx=5, pady=5)
    face_frame.grid(row=1, column=0)

    face = tk.Canvas(master=face_frame, bg="white", height=600, width=600)
    face_outline = face.create_oval(50, 50, 550, 550, fill="black")
    face.pack()
    
    plate_window.mainloop()

window_cfg = tk.Tk()
window_cfg.iconphoto(True, tk.PhotoImage(file="img/favicon.png"))
window_cfg.title("Octopus injector face pattern analysis")
window_cfg.columnconfigure(0, weight=1)
window_cfg.rowconfigure([0, 1], weight=1)
window_cfg.rowconfigure([2, 3], weight=0)

tab_parent = ttk.Notebook(window_cfg)
tab_manage = ttk.Frame(tab_parent)
tab_orifices = ttk.Frame(tab_parent)

tab_parent.add(tab_manage, text="Manage plates")
tab_parent.add(tab_orifices, text="Add orifices")
tab_parent.grid(row=0, column=0, sticky="w")

frame_top = tk.Frame(master=tab_manage, width=50, height=5, padx=5, pady=5)
frame_top.grid(row=1, column=0)

title = tk.Label(master=frame_top, text="Octopus injector face pattern analysis", height=1, width=50)
title.config(font=(96))
title.pack()

frame_diameter = tk.Frame(master=tab_manage, width=50, height=50, padx=5, pady=5)
frame_diameter.grid(row=2, column=0)

diameter_entry = tk.Entry(master=frame_diameter, width=10)
diameter_entry.grid(row=0, column=0, sticky="w")

diameter_entry_label = tk.Label(master=frame_diameter, text="plate diameter, mm")
diameter_entry_label.grid(row=0, column=1, sticky="e")

new_plate_button = tk.Button(master=tab_manage, text="Create injector plate", command=new_plate, width=30, height=3)
new_plate_button.grid(row=3, column=0)

window_cfg.mainloop()