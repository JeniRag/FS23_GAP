# -*- coding: utf-8 -*-
"""
Created on Fri Jun  2 15:22:07 2023

@author: Sujeni
"""

import numpy as np
import ase.units
import ase.io as ase_io
import os
import time

from sklearn.model_selection import KFold

from copy import deepcopy

from rascal.neighbourlist.structure_manager import (
                mask_center_atoms_by_species, mask_center_atoms_by_id)

from rascal.utils import from_dict, to_dict, CURFilter, dump_obj, load_obj, get_score, print_score, FPSFilter
from rascal.models import Kernel, sparse_points, train_gap_model, KRR, compute_KNM

from rascal.representations import SphericalExpansion, SphericalInvariants
from rascal.utils import (get_radial_basis_covariance, get_radial_basis_pca, 
                                  get_radial_basis_projections, get_optimal_radial_basis_hypers )

from rascal.utils import radial_basis
from utils import get_rmse


def write_log(message, logFile = "./Progress.log"):
    """
    To write logs with time points during computations.

    Parameters
    ----------
    message : str
        Message to write
    logFile : str, optional
        Files to write in it. The default is "./Progress.log".

    Returns
    -------
    None.

    """
    t=time.localtime()
    current_time = time.strftime("%H:%M:%S", t)
    f=open(logFile, "a")
    f.write(f"{current_time} {message} \n")
    f.close()
    

def reduce_forces(forces, idx):
    """
    Parameters
    ----------
    forces : np.ndarray
        3D array of forces with shape (N, M, 3).
    idx : np.ndarray
        indexes of the selected force components, with shape (N, Nforces)
        

    Returns
    -------
    forces_red :np.ndarray
        reduced forces with shape (N, Nforces, 3)

    """
    N = forces.shape[0]
    assert N==idx.shape[0]

    forces_reduced= forces[np.arange(N)[:, None], idx, :]     
    
    return forces_reduced



def extract_KNM(frames, kNM, index, Nforces):
    """
    Extract kNM values which belong to the corresponding frames in the given indeces.

    Parameters
    ----------
    frames : TYPE
        DESCRIPTION.
    kNM : np.ndarray
        Kernel Matrix of shape (N + N * Nforces * 3, M).
    index :np.ndarray
        1D Frame indexes to extract. 
    Nforces :int
        Number of Forces kept during building of model

    Returns
    -------
    np.ndarray
        Extracted kernel matrix

    """
    
    N= len(frames)
    assert kNM.shape[0] == N + N * Nforces * 3
    
    kNM_new_wog = kNM[:N][index]
    kNM_new_grad = np.concatenate([kNM[N + idx * Nforces * 3: N + (idx + 1) * Nforces * 3] for idx in index], axis=0)
    
    

    return np.concatenate([kNM_new_wog, kNM_new_grad], axis=0)

def ReduceFull_KNM(kNM, frames, force_idx):
    """
    Reduce Full KNM (with all force components) to one with reduced force components.

    Parameters
    ----------
    kNM : np.ndarray
        Kernel Matrix containing the similarity and gradients.
    frames : TYPE
        DESCRIPTION.
    force_idx : np.ndarray
        array of size (N, Nforces),

    Returns
    -------
    kNM_tot : TYPE
        DESCRIPTION.

    """
    
    # atoms per frame
    Natoms = len(frames[0])
    
    #number of forces kept
    Nforces = force_idx.shape[1]
    
    #Number of frames
    N = len(frames)

    kNM_new_grad = []
    
    kNM_wog = kNM[:N] #without gradient
    kNM_grad= kNM[N:] # with gradient

    for i in range(force_idx.shape[0]):
        start = i*Natoms*3 #skipy x,y,z of all the previous atoms
        for j in range(Nforces):
            s= start +force_idx[i,j] * 3
            e = s + 3 #includes x,y,z
            kNM_new_grad.append(kNM_grad[s: e])
            
    #kNM_new_grad = [kNM_grad[i * Natoms * 3 + force_idx[i, j] * 3:i * Natoms * 3 + force_idx[i, j] * 3 + 3] for i in range(force_idx.shape[0]) for j in range(Nforces)]

    
         
    kNM_new_grad = np.array(kNM_new_grad) #convert list to array
    kNM_new_grad = np.concatenate(np.array(kNM_new_grad), axis=0) #bring in right shape
    kNM_tot = np.concatenate([kNM_wog, kNM_new_grad], axis=0) #stack similarity values and new gradient together
    
    
    assert (kNM_tot.shape[0] == (N+ N*Nforces*3)) #check if shape is right
    return kNM_tot

def random_forces_index(frames, Nforces, seed=0):
    """
    Create random force indices.

    Parameters
    ----------
    frames : TYPE
        DESCRIPTION.
    Nforces : int
        Number of forcecomponents to keep.
    seed : int
        DESCRIPTION.

    Returns
    -------
    forces_idx : np.ndarray
        force indexes of shape (len(frames), Nforces)

    """
    np.random.seed(seed)
    N = len(frames)
   
    forces_idx = np.random.randint(0, len(frames[0]), size=(N, Nforces), dtype=int )
    
    return forces_idx


def predict_forces(model, frames, kNM, forces_index):
    """
    

    Parameters
    ----------
    model : TYPE
        GAP model.
    frames : TYPE
        DESCRIPTION.
    kNM : np.ndarray
        DESCRIPTION.
    forces_index : np.ndarray
        DESCRIPTION.

    Returns
    -------
    F : np.ndarray
        The predicted forces of size (len(frames), forces_index.shape[1])

    """
    
    forces_predicted = -(kNM[len(frames):] @ model.weights).reshape(-1,3)
    
    if forces_index is None: #when you use full force components
        forces_predicted= forces_predicted.reshape(-1, len(frames[0]), 3)
        
    else: #when you use reduced forces
        forces_predicted = forces_predicted.reshape(-1, forces_index.shape[1] ,3)
      
    return forces_predicted


def rs_to_gcm3(rs):
    return 2.6960431 / rs**3

        

    
