# -*- coding: utf-8 -*-
"""
Created on Mon Apr 17 09:42:50 2023

@author: Sujeni
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import re

from ase.visualize.plot import plot_atoms

import ase.units
import ase.io as ase_io

def PlotComparision(FileObject_list,  title=None,xlim=None, ylimE=None, ylimP=None, ylimT=None, T=8000):
    f, axs = plt.subplots(3, 1, sharex=True, figsize=(12, 6))
    f.subplots_adjust(hspace=0)
    if title !=None:
        f.suptitle(title)
   
    for F in FileObject_list:
        axs[0].plot(F.time, F.Energy, label=F.label)
        axs[1].plot(F.time, F.Pressure_tot, label=F.label)
        
        if np.any(type(F.Temperature) == np.ndarray):

            axs[2].plot(F.time, F.Temperature, label=F.label)
            axs[2].axhline(y=T, linestyle='--', c="black")
            
     
        
    axs[0].set_ylabel("Energy [EV]")
    axs[1].set_ylabel("Pressure [GPa]")
    axs[2].set_ylabel("Temperature [K]")  

    axs[2].set_xlabel("time [ps]")
    
    if xlim!=None:
        axs[1].set_xlim(xlim)
        axs[0].set_xlim(xlim)

    if ylimE!=None:
        axs[0].set_ylim(ylimE)
    if ylimP!=None:
        axs[1].set_ylim(ylimP)
    if ylimT!=None:
        axs[2].set_ylim(ylimT)
        
  
    axs[0].legend()
    axs[1].legend()
    axs[2].legend()
    
    return

def Plot(F, title=None, xlim=None, ylimT=None, ylimE=None, ylimP=None):
    f, axs = plt.subplots(4, 1, sharex=True, figsize=(12, 6))
    f.subplots_adjust(hspace=0)
    if title !=None:
        f.suptitle(title)
   

    axs[0].plot(F.time, F.Energy)
    axs[1].plot(F.time, F.Pressure_tot )

    if np.any(F.Temperature != None):
        axs[2].plot(F.time, F.Temperature)
        
    if F.file_type=="LAMMPS":
        #axs[3].plot(F.time,(F.PE- F.PE[0]), linewidth = 1.5, label = 'Change in Potential energy' )
        #axs[3].plot(F.time,(F.KE- F.KE[0]), linewidth = 1.5, label = 'Change in Kinetic energy' )
            
        axs[3].plot(F.time,(F.PE), linewidth = 1.5, label = 'Potential energy' )
        axs[3].plot(F.time,(F.KE), linewidth = 1.5, label = 'Kinetic energy' )
        
    axs[0].set_ylabel("Energy [EV]")
    axs[1].set_ylabel("Pressure [GPa]")
    axs[2].set_ylabel("Temperature [K]")  
    axs[3].set_ylabel("Energy [EV]")

    axs[3].set_xlabel("time [ps]")
    
    if xlim!=None:
        axs[1].set_xlim(xlim)
        axs[0].set_xlim(xlim)

    if ylimE!=None:
        axs[0].set_ylim(ylimE)
    if ylimP!=None:
        axs[1].set_ylim(ylimP)
    if ylimT!=None:
        axs[2].set_ylim(ylimT)
        
    axs[3].legend()
    
    return


def compare_atoms(file):
    f=open(file)
    line=f.readlines()
    match = re.search(r'CELL\(abcABC\):\s+([\d.]+)',line[1])
    # Check if a match was found
    if match:
        # Extract the float value
        cell_value = float(match.group(1))    
    f.close()

    plt.figure(figsize=(10,6))

    init = ase_io.read(file, index=0) #initial configuration
    plt.subplot(1,3,1)
    plt.title("initial configuration")
    plot_atoms(init)
    
    last=ase_io.read(file, index=-1) #last configuration
    plt.subplot(1,3,2)
    plt.title("last configuration")
    plot_atoms(last)
    
    plt.subplot(1,3,3)
    last.set_cell([cell_value, cell_value, cell_value, 90.0, 90.0, 90.0])
    last.set_pbc([True, True, True])
    last.wrap()
    plt.title("with periodic plot")
    plot_atoms(last)
    
