import matplotlib
matplotlib.use("TkAgg")

from tkinter import scrolledtext

from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)

from matplotlib.backends.backend_tkagg import \
    NavigationToolbar2Tk as NavigationToolbar2TkAgg

from tkinter.filedialog import askdirectory

from matplotlib.figure import Figure
import matplotlib.animation as animation
from matplotlib import style

import tkinter as tk
from tkinter import ttk

import urllib
import json

import pandas as pd
import numpy as np

from matplotlib import pyplot as plt
from gw_functions import *

LARGE_FONT= ("Verdana", 14)
NORM_FONT= ("Verdana", 10)
SMALL_FONT= ("Verdana", 8)
plt.style.use('ggplot')
import os

plt.rcParams['axes.facecolor'] = 'w'
plt.rcParams['grid.color'] = 'grey'
plt.rcParams['grid.alpha'] = 0.0
plt.rcParams['axes.linewidth'] = 0.5
plt.rc('axes', edgecolor='grey')

plt.rcParams['axes.spines.top'] = 0
plt.rcParams['axes.spines.right'] = 0
plt.rcParams['axes.spines.left'] = 1
plt.rcParams['axes.spines.bottom'] = 1
plt.rc('axes', edgecolor='grey')
plt.rcParams['image.cmap'] = 'Blues'

from itertools import cycle
lines = ["-", "--", "-.", ":"]
linecycler = cycle(lines)


f = Figure()
ax = f.add_subplot(111)


def popupmsg(msg):
    popup = tk.Tk()
    popup.wm_title("!")
    label = ttk.Label(popup, text=msg, font=LARGE_FONT)
    label.pack(side="top", fill="x", pady=10)
    B1 = ttk.Button(popup, text="Okay", command = popup.destroy)
    B1.pack()

    popup.mainloop()

     

class GW_gui(tk.Tk):

    def __init__(self, *args, **kwargs):
        
        tk.Tk.__init__(self, *args, **kwargs)

        tk.Tk.iconbitmap(self)
        tk.Tk.wm_title(self, "Giraldez and Woolhiser")
        
        container = tk.Frame(self)
        container.pack(side="top", fill="both", expand = True)
        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)


        # menubar = tk.Menu(container)
        # filemenu = tk.Menu(menubar, tearoff=0)
        # filemenu.add_command(label="Save settings", 
        #     command = lambda: popupmsg("Not supported just yet!"))
        # filemenu.add_separator()
        # filemenu.add_command(label="Exit", command=quit)
        # menubar.add_cascade(label="File", menu=filemenu)

        # tk.Tk.config(self, menu=menubar)

        frame = graph_page(container, self)
        frame.grid(row=0, column=0, sticky="nsew")

    
def test_valid(param):
    if param["ksatV"] <= 0:
        return 0
    else:
        return True

class graph_page(tk.Frame):

    def reset(self):

        self.action.configure(text='Run the model' )

        self.clear.configure(state = tk.DISABLED)  
        self.save.configure(state = tk.DISABLED)  

        self.scr.delete('1.0', tk.END)

        ax.clear()

        self.canvas.draw()

    def clickMe(self):
            
        Ks = float(self.Ks.get())
        if self.Ks_unit.get() == "cm/hr":
            Ks = Ks/3.6e5

        Ao = float(self.Ao.get())
        if self.Ao_unit.get() == "cm/sqrt(hr)":
            Ao = Ao/100./60.
        
        rain = float(self.rain.get())
        if self.rain_unit.get() == "cm/hr":
            rain = rain/3.6e5

        param = {'ksatV': Ks,
             'rain': rain,
             'So': 0.1,
             'eta' : 0.5,
             'alpha' : 0.1,
             'a': 2./3,
             't_rain': 1200,
             'L': 50.0,
             'dt': 1.0,
             'ntstep': 500,
             'delta_theta': 0.24,
             'Ao': Ao}
 
        self.valid_param = test_valid(param)

        value = "Ks = {0} {1} \n".format(self.Ks.get(), self.Ks_unit.get())
        value += "Ao = {0} {1} \n".format(self.Ao.get(), self.Ao_unit.get())             
        if not self.valid_param:
            popupmsg("Double check those parameters!")
            self.scr.insert(tk.INSERT, "Invalid parameters:" + '\n')   
            self.scr.insert(tk.INSERT, value + '\n')   

        if self.valid_param:
 
            self.res = Comparison_function2(param)
            self.res = pd.Series(self.res)
            
            self.count += 1
            self.scr.insert(tk.INSERT, 
                "Model run {0} parameters:".format(self.count) + '\n')   
            self.scr.insert(tk.INSERT, value + '\n')   
            
            self.action.configure(text='Run again!' )
            self.plot(ax)

        self.clear.configure(state = tk.NORMAL)
        self.save.configure(state = tk.NORMAL)  

    def plot(self,ax):
        
        x =  np.arange(200)
        y = np.random.randint(0, 200, 200)
        # ax.clear()

        ax.plot(self.res.t/60, self.res.q/self.res.L*3.6e5, label = self.count)

        ax.legend(bbox_to_anchor=(0.9, 0.8, 1, .102), loc=3,
                 ncol=1, borderaxespad=0)

        ax.set_title("Hydrograph")    
        ax.set_xlabel("minutes")

        self.canvas.draw()

    def save_file(self):

        tk.Tk().withdraw() 
        filepath = askdirectory() 
   
        filepath = os.path.join(filepath, "GW_output.xlsx")

        res = self.res
        univariate = ['ksatV', 'rain', 'So', 'eta', 'alpha', 'a', 't_rain', 
                    'L', 'ntstep', 'delta_theta', 'Ao', 'Kr', 't_pond']  
        
        df = pd.DataFrame({"t" :res["t"], "q" : res["q"]})

        with pd.ExcelWriter(filepath) as writer:  
            res[univariate].to_excel(writer, sheet_name='univariate')
            df.to_excel(writer, sheet_name='hydrograph')      
        # self.res.to_csv(filepath)


    def __init__(self, parent, controller):

        
        self.valid_param = 0
        self.count = 0

        tk.Frame.__init__(self, parent)
        label = tk.Label(self, text="The Giraldez and Woolhiser GUI!", 
            font=LARGE_FONT)
        label.pack(pady=10,padx=10)
       
        # Create a container to hold labels
        labelsFrame = ttk.LabelFrame(self, text=' Select model parameters ')
        labelsFrame.pack(expand=True) 

        rain_row = 0
        Ks_row = 1
        Ao_row = 2
        # Place labels into the container element - vertically
        ttk.Label(labelsFrame, text="rain = ").grid(column=0, row=0)

        self.rain = tk.StringVar()
        choose_rain = ttk.Entry(labelsFrame, width=12, 
            textvariable = self.rain)
        self.rain.set(5.0)
        choose_rain.grid(column=1, row=rain_row, sticky='W')

        self.rain_unit = tk.StringVar()
        choose_rain_unit = ttk.Combobox(labelsFrame, width=12, 
            textvariable = self.rain_unit)

        choose_rain_unit['values'] = ("cm/hr", "m/s")
        choose_rain_unit.grid(column=2, row=rain_row)
        choose_rain_unit.current(0)

        ttk.Label(labelsFrame, text="Ks = ").grid(column=0, row= Ks_row)

        self.Ks = tk.StringVar()
        choose_Ks = ttk.Entry(labelsFrame, width=12, 
            textvariable = self.Ks)
        self.Ks.set(1.0)
        choose_Ks.grid(column=1, row=Ks_row, sticky='W')

        self.Ks_unit = tk.StringVar()
        choose_Ks_unit = ttk.Combobox(labelsFrame, width=12, 
            textvariable = self.Ks_unit)

        choose_Ks_unit['values'] = ("cm/hr", "m/s")
        choose_Ks_unit.grid(column=2, row=Ks_row)
        choose_Ks_unit.current(0)

        ttk.Label(labelsFrame, text="Ao = ").grid(column=0, row=Ao_row)

        self.Ao = tk.StringVar()
        self.Ao.set(1.5)

        choose_Ao = ttk.Entry(labelsFrame, width=12, 
            textvariable=self.Ao)
        choose_Ao.grid(column=1, row=Ao_row, sticky='W')


        self.Ao_unit = tk.StringVar()
        choose_Ao_unit = ttk.Combobox(labelsFrame, width=12, 
            textvariable = self.Ao_unit)
        choose_Ao_unit['values'] = ("cm/sqrt(hr)", "m/sqrt(s)")
        choose_Ao_unit.grid(column=2, row=Ao_row)
        choose_Ao_unit.current(0)


        self.action = ttk.Button(labelsFrame, 
            text="Run the model",  command=self.clickMe)
        self.action.grid(column=0, row=6)   
        
        self.clear = ttk.Button(labelsFrame, 
            text="Clear results",  command=self.reset, 
            state=tk.DISABLED)
        self.clear.grid(column=1, row=6)  

        self.save = ttk.Button(labelsFrame, 
            text="Save results",  
            command=self.save_file, 
            state=tk.DISABLED)
        self.save.grid(column=2, row=6)  

        scrolW  = 30; scrolH  =  10
        self.scr = scrolledtext.ScrolledText(labelsFrame, 
            width=scrolW, height=scrolH, wrap=tk.WORD)
        self.scr.grid(column=0, row=10, sticky='WE', columnspan=3)



        # Add some space around each label
        for child in labelsFrame.winfo_children():
            child.grid_configure(padx=10, pady = 4)


        self.canvas = FigureCanvasTkAgg(f, self)
        self.canvas.get_tk_widget().pack(side=tk.BOTTOM, 
            fill=tk.BOTH, expand=True)
        self.canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        toolbar = NavigationToolbar2TkAgg(self.canvas, self)
        toolbar.update()

app = GW_gui()
app.geometry("1000x800")
app.mainloop()