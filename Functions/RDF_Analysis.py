# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 11:52:34 2023

@author: Sujeni
"""

from tqdm.notebook import tnrange
import numpy as np
import os
import ase.io as ase_io
from ase.visualize import view
from ase.visualize.plot import plot_atoms
from ase.geometry import get_distances
from ase.geometry.rdf import get_rdf
from ase.geometry.analysis import Analysis
import matplotlib.pyplot as plt

import re

def set_periodic(file_path, index=1):
    
    """
    file_path: string
        Path of the file to read.
    index: int, slice or str
        slicing indexes for the fames
    
    Returns ase cell, with coordinates set to periodic. We assume here that the cell sizes of frames belonging to the same file are constant.
    
    """
    
    #read the frames
    frames=ase_io.read(file_path, index=index)
    
    num_lines = 10
    with open(file_path, "r") as file:
        line = []
        for _ in range(num_lines):
            l = file.readline()
            if not l:
                break  # Break the loop if end of file is reached
            line.append(l)
    
    
    cell_value = None
    match = re.search(r'CELL\(abcABC\):\s+([\d.]+)', line[1]) #search for cell value
    
    # Check if a match was found
    if match:
        # Extract the float value
        cell_value = float(match.group(1))
        
    else:
        cell_value = np.round(frames[0].cell.lengths()[0],5)
    
    assert(cell_value != 0.0)
    if cell_value == None:
        raise ValueError("cell dimension not detected")

    for i in tnrange(len(frames)):
        fr = frames[i]
        fr.set_cell([cell_value, cell_value, cell_value, 90.0, 90.0, 90.0])
        fr.set_pbc([True, True, True])
        fr.wrap()
    
    return frames
    
    
def RDF(file_path, index=1, nbins=100, rmax=None):
    """
        file_path: str
            path to the tile
        index: int, slice or str
            slicing index of the frames
        nbins:int
            number of nbins
        rmax: float
            maximum possible radius for rdf calculation
        save_path: str
            path to store rdf numpy array
        
        
        returns numpy array containing [r, rdf]
    """

 

    #make frames periodic
    frames = set_periodic(file_path, index)
    
    #empty arrays to store values
    gr_list = np.zeros((len(frames), nbins))
    r_list=np.zeros((len(frames), nbins))
    
    #determine maximum possible rmax
    if rmax==None:
        rmax = frames[0].cell.lengths()[0]/2*0.98
        
    for i in tnrange(len(frames)):
        atoms = frames[i]
        
        # Compute the radial distribution function
        g_r, r = get_rdf(atoms, rmax=rmax, nbins=nbins)
        gr_list[i,:]=g_r
        r_list[i,:]=r #will be the same value

    #calculate the average
    g_r_mean=np.mean(np.array(gr_list),axis=0)
    
    #plot
    # plt.plot(r_list[0,:], g_r_mean)
    # plt.title("RDF")
    # plt.xlabel("r [Angstrom]")
    # plt.ylabel("g(r)")

        
    return  [r_list[0,:], g_r_mean]
    