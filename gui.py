from octopus import Fluid, Orifice, PropertySource, Manifold, Element
import tkinter as tk

window_cfg = tk.Tk()
window_cfg.title("Octopus - injector face pattern analysis")
frame_top = tk.Frame(master=window_cfg, width=200, height=5, padx=5, pady=5)
frame_top.grid(row=0, column=0)

title = tk.Label(master=frame_top, text="Octopus - injector face pattern analysis", height=1)
title.pack()

window_cfg.mainloop()