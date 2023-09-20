# -*- coding: utf-8 -*-
"""
Created on Mon Apr 17 12:42:02 2023

@author: Sujeni
"""
import ase.units

import ase.io as ase_io
import re

import os

def write_periodic(file_path):
    """Write outputs of periodic boundary condition simulations, 
    performed with Lammps, to periodice coordinates."""
    
    frames=ase_io.read(file_path, index=":")
    
    f=open(file_path)
    line=f.readlines()
    match = re.search(r'CELL\(abcABC\):\s+([\d.]+)',line[1])
    
    # Check if a match was found
    if match:
        # Extract the float value
        cell_value = float(match.group(1))
    else:
        cell_value = frames[0].cell.lengths()[0]
        
    assert((cell_value != 0.0 or cell_value != None ))
    
    f.close()

    for fr in frames:
        fr.set_cell([cell_value, cell_value, cell_value, 90.0, 90.0, 90.0])
        fr.set_pbc([True, True, True])
        fr.wrap()
    
    head, tail = os.path.split(file_path)
    ase_io.write(head+"/periodic.xyz", frames, format="extxyz")
        
    return

def read_iPi(file_path, index = ":"):
    
    frames=ase_io.read(file_path, index=index)
    #f=open(file_path, "r")
    #lines=f.readlines()
    #f.close()
    
  
    num_lines = 10

    with open(file_path, "r") as file:
        lines = []
        for _ in range(num_lines):
            line = file.readline()
            if not line:
                break  # Break the loop if end of file is reached
            lines.append(line)
    
    #print(lines)
    
    cell_value = None
    match = re.search(r'CELL\(abcABC\):\s+([\d.]+)', lines[1]) #assume volume is constant
    
    # Check if a match was found
    if match:
        # Extract the float value
        cell_value = float(match.group(1))
        
    else:
        cell_value = frames[0].cell.lengths()[0]
    
    assert(cell_value != 0.0)
    if cell_value == None:
        raise ValueError("cell dimension not detected")
        
    for i in range(len(frames)):
        fr = frames[i]
        fr.set_cell([cell_value, cell_value, cell_value, 90.0, 90.0, 90.0])    
        
    return frames

        

def rewrite_files_for_dynamics(file, destination_dir):
    """
    Takes an output file of QE simulation, extracts the first coordinates, and saves it in an format, that is used in NVT simulations using i-Pi.

    Parameters
    ----------
    file :str
        File path pointing to the output file of the DFT-reference simulation.
    destination_dir : str
        File path to store the i-Pi input files.

    Returns
    -------
    None.

    """
    frames=ase_io.read(file, index="0") #converts it arleady into armstrong
    L=frames.cell.lengths()[0] #converts in armstrong
    
    if not os.path.exists(destination_dir):
        os.makedirs(destination_dir)

    angle=90.00000
    
    Nframes = len(frames)    
    f = open(destination_dir +f"pos.xyz", "w")
    f.write(f"{len(frames)} \n")
    f.write(f"# CELL(abcABC): {L:.5f} {L:.5f} {L:.5f} {angle:.5f} {angle:.5f} {angle:.5f} Step: 0 Bead: 0 positions" + "{angstrom} cell{angstrom}\n")
    for i in range(Nframes):
        f.write(f"{frames.get_chemical_symbols()[i]} {frames.positions[i][0]:.5e} {frames.positions[i][1]:.5e} {frames.positions[i][2]:.5e}\n")
    f.close()
    return