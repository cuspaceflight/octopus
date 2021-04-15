from octopus import Fluid, Orifice, PropertySource, Manifold, Element
import tkinter as tk

window_cfg = tk.Tk()
window_cfg.iconphoto(True, tk.PhotoImage(file="img/favicon.png"))
window_cfg.title("Octopus - injector face pattern analysis")
window_cfg.columnconfigure(0, weight=1)
window_cfg.rowconfigure([0, 1], weight=1)
frame_top = tk.Frame(master=window_cfg, width=200, height=5, padx=5, pady=5)
frame_top.grid(row=0, column=0)

title = tk.Label(master=frame_top, text="Octopus - injector face pattern analysis", height=1)
title.pack()

frame_bottom = tk.Frame(master=window_cfg, height=50, padx=5, pady=5)
frame_bottom.grid(row=1, column=0)

diameter_entry = tk.Entry(master=frame_bottom, width=10)
diameter_entry.grid(row=0, column=1, sticky="w")

diameter_entry_label = tk.Label(master=frame_bottom, text="plate diameter, mm")
diameter_entry.grid(row=0, column=1, sticky="e")
window_cfg.mainloop()