# -*- coding: utf-8 -*-
"""
Created on Wed May 17 11:04:28 2023

@author: Sujeni
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pandas as pd
import os
import re


def read_dt_ref(path):
    """
    read time step from file and return dictionary
    """
   
    dt_ref = {} #dt of reference
    files = os.listdir(path)
    
    for file in files:
        if file.endswith(".in"): #files ending with .in
            pattern =  r"rho(\d+\.\d+)"
            match = re.search(pattern, file)
            if match:
                rho = float(match.group(1)) #extract rho of file name
                
            f = open(path+file, "r")
            lines = f.readlines()
            pattern = r'^\s*dt\s*=\s*([\d\.]+)d0' #extract time step [a.u]
            
            for line in lines:
                match = re.match(pattern, line)
                if match:
                    dt_ref[rho] = float(match.group(1)) #save in dictionary
                    break
        
            f.close()

    return dt_ref


def InspectionPlot(F, F_ref=None, xlim = None, xlim_ref = None, ylimP=None):
    
    if F_ref !=None:
        fig, axs = plt.subplots(2, 1, figsize=(16, 4), sharex=True)
        fig.subplots_adjust(hspace=0)
        fig.suptitle("Reference")
        
        #axs[0].set_title("Energy")
        axs[0].set_ylabel("Energy [eV]")
        axs[0].set_xlabel("time [ps]")
        axs[0].plot(F_ref.time, F_ref.Energy)
        
    
        #axs[1].set_title("Pressure vs. Time of Reference")
        axs[1].set_ylabel("Pressure [GPa]")
        axs[1].set_xlabel("time [ps]")
        axs[1].plot(F_ref.time, F_ref.Pressure_tot)
        if ylimP is not None:
            axs[1].set_ylim(ylimP)
        
        if xlim_ref!=None:
            for i,values in enumerate([F_ref.Energy, F_ref.Pressure_tot]):
                axs[i].vlines(xlim_ref[0], values.min(), values.max(), alpha = 0.3)
                axs[i].vlines(xlim_ref[1], values.min(), values.max(), alpha = 0.3)
        
        fig, axs = plt.subplots(2, 1, figsize=(16, 4), sharex=True)
        fig.subplots_adjust(hspace=0)
        fig.suptitle("Reference Zoomed In")
        
        #axs[0].set_title("Energy")
        axs[0].set_ylabel("Energy [eV]")
        axs[0].set_xlabel("time [ps]")
        axs[0].plot(F_ref.time, F_ref.Energy)
        
        #axs[1].set_title("Pressure vs. Time of Reference")
        axs[1].set_ylabel("Pressure [GPa]")
        axs[1].set_xlabel("time [ps]")
        axs[1].plot(F_ref.time, F_ref.Pressure_tot)
    
        if xlim_ref!=None:
            axs[1].set_xlim(xlim_ref)
            axs[0].set_xlim(xlim_ref)
            stepsize=(xlim_ref[1]-xlim_ref[0])/20
            axs[1].set_xticks(np.arange(xlim_ref[0], xlim_ref[1], stepsize))
            
        
        
    fig, axs = plt.subplots(3, 1, figsize=(16, 6), sharex=True)
    fig.subplots_adjust(hspace=0)
    plt.tight_layout()
    #plt.autoscale(enable=True, axis='y', tight=False)
    
    fig.suptitle("GAP")
    
    i=0
    #axs[i].set_title("Energy vs. Time Lammps")
    axs[i].set_ylabel("Energy [eV]")
    #axs[i].set_xlabel("time [ps]")
    axs[i].plot(F.time, F.Energy)
    
    i=1
    #axs[i].set_title("Pressure vs. Time Lammps")
    axs[i].set_ylabel("Pressure [GPa]")
    #axs[i].set_xlabel("time [ps]")
    axs[i].plot(F.time, F.Pressure_tot)
    
    i=2
    #axs[i].set_title("Temperature vs. Time of Lammps")
    axs[i].set_ylabel("Temperature [ps]")
    axs[i].set_xlabel("time [ps]")
    axs[i].plot(F.time, F.Temperature)   
    
    if xlim!=None:
        for i,values in enumerate([F.Energy, F.Pressure_tot, F.Temperature]):
            axs[i].vlines(xlim[0], values.min(), values.max(), alpha = 0.5)
            axs[i].vlines(xlim[1], values.min(), values.max(), alpha = 0.5)
    plt.show()
       
        

    fig, axs = plt.subplots(3, 1, figsize=(16, 6), sharex=True)
    fig.subplots_adjust(hspace=0)
    plt.tight_layout()
    #plt.autoscale(enable=True, axis='y', tight=False)
    
    fig.suptitle("GAP Zoomed in")
    
    i=0

    axs[i].set_ylabel("Energy [eV]")
    axs[i].plot(F.time, F.Energy)
    
    i=1
    axs[i].set_ylabel("Pressure [GPa]")
    axs[i].plot(F.time, F.Pressure_tot)
    
    i=2
    axs[i].set_ylabel("Temperature [ps]")
    axs[i].set_xlabel("time [ps]")
    axs[i].plot(F.time, F.Temperature)

    if xlim!=None:
        axs[0].set_xlim(xlim)
        axs[1].set_xlim(xlim)
        axs[2].set_xlim(xlim)
        stepsize=(xlim[1]-xlim[0])/20
        axs[2].set_xticks(np.arange(xlim[0], xlim[1], stepsize))
        mask =  ((F.time >= xlim[0]) & (F.time < xlim[1]))
        
        axs[0].set_ylim(min(F.Energy[mask]), max(F.Energy[mask]) )
        axs[1].set_ylim(min(F.Pressure_tot[mask]), max(F.Pressure_tot[mask]) )
        axs[2].set_ylim(min(F.Temperature[mask]), max(F.Temperature[mask]))
        
def standard_error(array):
    return np.std(array, ddof=1) / np.sqrt(len(array))    
    
def decorrelated(LN, RN=None, F=None, F_ref=None, start =0, end = None, n_atoms= None, start_ref=0, plot = True):
    """
    

    Parameters
    ----------
    LN :int
        Points of GAP simulation
    RN : int, optional
        Points of reference. The default is None.
    F : FileObject, optional
        DESCRIPTION. The default is None.
    F_ref : FileObject, optional
        DESCRIPTION. The default is None.
    start : int, optional
       starting point. The default is 0.
    end : int, optional
       ending point. The default is None.
    n_atoms : int, optional
        number of atoms. The default is None.
    start_ref : TYPE, optional
        DESCRIPTION. The default is 0.
    plot : TYPE, optional
        DESCRIPTION. The default is True.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
   
    eV_to_Hartree = 0.0367
    
    #filter desired points
    LEnergy = F.Energy[start:end:LN] * eV_to_Hartree
    LPressure = F.Pressure_tot[start:end:LN]
    LTemperature=F.Temperature[start:end:LN]
    N = len(LPressure)
    print(f"length decimated array Lammps: {N}")
    
    #calculate mean and std
    mean_LPressure = np.mean(LPressure)
    sigma_LPressure = standard_error(LPressure)
    
    mean_LEnergy = np.mean(LEnergy)
    sigma_LEnergy = standard_error(LEnergy)
    
    mean_LTemperature = np.mean(LTemperature)
    sigma_LTemperature = standard_error(LTemperature)
    
    #if reference contained
    if F_ref != None:    
    
        REnergy = F_ref.Energy[start_ref::RN] * eV_to_Hartree
        RPressure=F_ref.Pressure_tot[start_ref::RN]
        print(f"length decimated array Ref: {len(RPressure)}")
    
        mean_RPressure = np.mean(RPressure)
        sigma_RPressure = standard_error(RPressure)
            
        mean_REnergy = np.mean(REnergy)
        sigma_REnergy = standard_error(REnergy)
        
        #print information on screen
        print(f" Pressure gap: {mean_LPressure:.2f} +/- { sigma_LPressure:.2f} GPa")
        print(f" Pressure ref: {mean_RPressure:.2f}  +/-  {sigma_RPressure:.2f} GPa")
               
        print(f"Temperature gap: {np.mean(LTemperature):.2f} +/- {sigma_LTemperature:.2f} K")
        
        print(f" Energy/atom gap: {mean_LEnergy:.2f} +/- {sigma_LEnergy:.2f} Hartree")
        print(f" Energy/atom ref: {mean_REnergy:.2f} +/- {sigma_REnergy:.2f} Hartree")
        
        if plot:
            fig, axs = plt.subplots(3,1,figsize=(16, 6), sharex=True)
            
            axs[0].set_title("Decorrelated Energy")
            axs[0].plot(LEnergy, label = "gap")
            axs[0].plot(REnergy, label = "ref")
            axs[0].set_ylabel("Energy [eV]")
            
            axs[1].set_title("Decorrelated Pressure")
            axs[1].plot(LPressure, label="gap")
            axs[1].plot( RPressure, label="ref")
            axs[1].set_ylabel("Pressure [GPa]")
            
            i=2
            axs[i].set_title("Decorrelated Temperature")
            axs[i].plot(LTemperature)
            axs[i].set_ylabel("Temperature [K]")
             
        return [mean_LPressure, sigma_LPressure], [mean_RPressure, sigma_RPressure], [mean_LEnergy, sigma_LEnergy], [mean_REnergy, sigma_REnergy], [mean_LTemperature, sigma_LTemperature], N

    else:

        print(f" Pressure gap: {mean_LPressure:.2f} +/- { sigma_LPressure:.2f} GPa")   
        print(f"Temperature gap: {np.mean(LTemperature):.2f} +/- {sigma_LTemperature:.2f} K")        
        print(f" Energy/atom gap: {mean_LEnergy:.2f} +/- {sigma_LEnergy:.2f} Hartree")
        
        if plot:
            fig, axs = plt.subplots(3,1,figsize=(16, 6), sharex=True)
            
            axs[0].set_title("Decorrelated Energy")
            axs[0].plot(LEnergy, label = "gap")
          
            axs[0].set_ylabel("Energy [eV]")
            
            axs[1].set_title("Decorrelated Pressure")
            axs[1].plot(LPressure, label="gap")
        
            axs[1].set_ylabel("Pressure [GPa]")
            
            i=2
            axs[i].set_title("Decorrelated Temperature")
            axs[i].plot(LTemperature)
            axs[i].set_ylabel("Temperature [K]")
            
            plt.legend()
              

        return [mean_LPressure, sigma_LPressure], [mean_LEnergy, sigma_LEnergy], [mean_LTemperature, sigma_LTemperature], N

            
def ShowSelectedPoints(F, N, start = 0, end = None, xlim=None):
    
    if type(F.Temperature)==np.ndarray:
        
        f, axs = plt.subplots(3, 1, sharex=True, figsize=(16, 6))
        f.subplots_adjust(hspace=0)
        f.suptitle(F.label )
        
        axs[0].plot(F.time, F.Energy, alpha = 0.5, c="blue")
        axs[0].scatter(F.time[start:end:N], F.Energy[start:end:N], c="red")
        axs[0].set_ylabel("Energy [eV]")
        
        axs[1].plot(F.time, F.Pressure_tot, alpha=0.5, c= "blue")
        axs[1].scatter(F.time[start:end:N], F.Pressure_tot[start:end:N], c="red")
        axs[1].set_ylabel("Pressure [GPa]")

        axs[2].plot(F.time, F.Temperature, alpha=0.5, c= "blue")
        axs[2].scatter(F.time[start:end:N], F.Temperature[start:end:N], c="red")
        axs[2].set_ylabel("Temperature [K]")
    
        axs[2].set_xlabel("time [ps]")
        

        
        if xlim is not None:
            for i, values in enumerate([F.Energy, F.Pressure_tot, F.Temperature]):
                axs[i].vlines(xlim[0], values.min(), values.max(), alpha = 0.5)
                axs[i].vlines(xlim[1], values.min(), values.max(), alpha = 0.5)
        
            f, axs = plt.subplots(3, 1, sharex=True, figsize=(16, 6))
            f.subplots_adjust(hspace=0)
            f.suptitle(F.label +" Zoom" )

            axs[0].plot(F.time, F.Energy, alpha = 0.5, c="blue")
            axs[0].scatter(F.time[start:end:N], F.Energy[start:end:N], c="red")
            axs[0].set_ylabel("Energy [eV]")

            axs[1].plot(F.time, F.Pressure_tot, alpha=0.5, c= "blue")
            axs[1].scatter(F.time[start:end:N], F.Pressure_tot[start:end:N], c="red")
            axs[1].set_ylabel("Pressure [GPa]")

            axs[2].plot(F.time, F.Temperature, alpha=0.5, c= "blue")
            axs[2].scatter(F.time[start:end:N], F.Temperature[start:end:N], c="red")
            axs[2].set_ylabel("Temperature [K]")

            axs[2].set_xlabel("time [ps]")


            axs[2].set_xlim(xlim)
            
            stepsize=(xlim[1]-xlim[0])/20
            axs[1].set_xticks(np.arange(xlim[0], xlim[1], stepsize))
            mask =  ((F.time >= xlim[0]) & (F.time < xlim[1]))
            
            axs[0].set_ylim(min(F.Energy[mask]), max(F.Energy[mask]) )
            axs[1].set_ylim(min(F.Pressure_tot[mask]), max(F.Pressure_tot[mask]) )
            axs[2].set_ylim(min(F.Temperature[mask]), max(F.Temperature[mask]))
            
            
        
    else:
        
        f, axs = plt.subplots(2, 1, sharex=True, figsize=(16, 4))
        f.subplots_adjust(hspace=0)
        plt.suptitle(F.label)
        
        axs[0].plot(F.time, F.Energy, alpha=0.5, c = "blue")
        axs[0].scatter(F.time[start:end:N], F.Energy[start:end:N], c="red")
        axs[0].set_ylabel("Energy [eV]")
        
        axs[1].plot(F.time, F.Pressure_tot, alpha=0.5, c= "blue")
        axs[1].scatter(F.time[start:end:N], F.Pressure_tot[start:end:N], c="red")
        axs[1].set_ylabel("Pressure [GPa]")
        
        axs[1].set_xlabel("time [ps]")
        
        if xlim != None:
            for i, values in enumerate([F.Energy, F.Pressure_tot]):
                axs[i].vlines(xlim[0], values.min(), values.max(), alpha = 0.5)
                axs[i].vlines(xlim[1], values.min(), values.max(), alpha = 0.5)
            
            f, axs = plt.subplots(2, 1, sharex=True, figsize=(16, 4))
            f.subplots_adjust(hspace=0)
            plt.suptitle(F.label+ " Zoom")

            axs[0].plot(F.time, F.Energy, alpha=0.5, c = "blue")
            axs[0].scatter(F.time[start:end:N], F.Energy[start:end:N], c="red")
            axs[0].set_ylabel("Energy [eV]")

            axs[1].plot(F.time, F.Pressure_tot, alpha=0.5, c= "blue")
            axs[1].scatter(F.time[start:end:N], F.Pressure_tot[start:end:N], c="red")
            axs[1].set_ylabel("Pressure [GPa]")

            axs[1].set_xlabel("time [ps]")
            
            
            axs[0].set_xlim(xlim)
            axs[1].set_xlim(xlim)
            
            
            

    print(f"{len(F.Pressure_tot[start:end:N])} elements")

def FindN(dt, lambd):
    #dt: time step of simulation in ps
    #lambd: time in ps

    N=lambd / (dt)
    print(N)
    
    
def plot_RDF(d_T, T, rows, cols, width = 16, height = 12, sup_gap = 0.95):
    rho_list = list(d_T.keys())

    
    fig, axs = plt.subplots(rows, cols, figsize = (width, height), sharex="col", sharey=True)
    
    if rows==1:
        axs = np.atleast_2d(axs)

    fig.suptitle(f"RDF for different densities at {T} K", y= sup_gap, size = 16)
    #fig.tight_layout()

    for idx, rho in enumerate(rho_list):
        i = idx // cols 
        j = idx % cols

        #print(f"{idx}: [{i}. {j}]")

        rdf_path=f"../dynamics/outputs/{model}/T{T}K/rho{rho}/rdf.npy"
        rdf = np.load(rdf_path)
        axs[i][j].plot(rdf[0,:], rdf[1,:], label="lammps", c="dodgerblue")

        rdf_path = f"../dynamics/inputs/T{T}K/rho{rho}/rdf.npy"
        rdf=np.load(rdf_path)
        axs[i][j].plot(rdf[0,:], rdf[1,:], label="ref", c="darkorange")




        axs[i][j].set_title(f" \n rho = {rho} $g/cm^3$", size = 14)



    for i in range(rows):
        for j in range(cols):
            axs[i][0].set_ylabel("RDF", size =13)
            axs[rows-1][j].set_xlabel("r [Angstrom]", size = 13)
        
       

    axs[0][0].legend()
    
def create_Table(d_T, T, round_value=1):
    data_list=[]
    rho_list = list(d_T.keys())
    for rho in rho_list:
        Lammps_P, Ref_P, Lammps_E, Ref_E, Lammps_T, N = d_T[rho]
        data = [Lammps_P[0], Lammps_P[1], Ref_P[0], Ref_P[1], Lammps_E[0], Lammps_E[1], Ref_E[0], Ref_E[1], Lammps_T[0], Lammps_T[1], N]
        data = [np.round(d,round_value) for d in data]
        data_list.append(data)
    
    df = pd.DataFrame(data_list, columns = ["LP", "sigma_LP", "RP", "sigma_RP","LE", "sigma_LE", "RE", "sigma_RE", "LT","sigma_LT", "N"], index=rho_list, dtype = float)

    df_new = pd.DataFrame()
    
    df_new['Pressure [GPa]'] = df["LP"].astype(str) + ' $\pm$ ' + df['sigma_LP'].astype(str)
    df_new['Pressure Ref [GPa]'] = df["RP"].astype(str) + ' $\pm$ ' + df['sigma_RP'].astype(str)
    df_new['Energy/atom [Hartrees]'] = df["LE"].astype(str) + ' $\pm$ ' + df['sigma_LE'].astype(str)
    df_new['Energy/atom Ref [Hartrees]'] = df["RE"].astype(str) + ' $\pm$ ' + df['sigma_RE'].astype(str)

    #df_new["Temperature [K]"] = df["LT"].astype(str) + " $\pm$ " + df["sigma_LT"].astype(str)
    df_new = df_new.rename_axis("density [g/cm3]", axis=0)
   
    return df, df_new


def create_Table_noRef(d_T, T, round_value=1):
    """
    

    Parameters
    ----------
    d_T : Dataframe
        Dataframe with density, Pressure, Energy values.
    T : int
        Temperature

    Returns
    -------
    df : dataframe
        DESCRIPTION.
    df_new : dataframe
        DESCRIPTION.

    """
    data_list=[]
    rho_list = list(d_T.keys())
    for rho in rho_list:
        Lammps_P, Lammps_E, Lammps_T, N = d_T[rho]
        data = [Lammps_P[0], Lammps_P[1], Lammps_E[0], Lammps_E[1],  Lammps_T[0], Lammps_T[1], N]
        data = [np.round(d,round_value) for d in data]
        data_list.append(data)
    
    df = pd.DataFrame(data_list, columns = ["LP", "sigma_LP", "LE", "sigma_LE", "LT","sigma_LT", "N"], index=rho_list, dtype = float)

    df_new = pd.DataFrame() 
    
    df_new['Pressure [GPa]'] = df["LP"].astype(str) + ' $\pm$ ' + df['sigma_LP'].astype(str)    
    df_new['Energy/atom [Hartrees]'] = df["LE"].astype(str) + ' $\pm$ ' + df['sigma_LE'].astype(str)
    df_new = df_new.rename_axis("density [g/cm3]", axis=0)
    return df, df_new

def plot_physics(df, df_noRef, T, zoomsize=None, xlimit=None, ylimit=None, reference = None):
    
    fig, axs = plt.subplots(1,2,figsize=(12, 5))
    if df_noRef is not None:
        df = pd.concat([df, df_noRef], sort = False)
    
    df=df.sort_index()
    x = list(df.index)
    
    dfRef = df[df['RE'].notna()]
    xRef = list(dfRef.index)
    

    lw = 0.5
    linestyle="--"
    
    #Gap
    axs[0].errorbar(x=x, y = df["LP"], yerr=df["sigma_LP"], marker="o", markersize=3, uplims=True, linestyle=linestyle, lolims=True, lw=lw, elinewidth=0.5, label ="GAP", capsize=0.5, c = "blue")
    #DFT-ref
    axs[0].errorbar(x=xRef, y = dfRef["RP"], yerr=dfRef["sigma_RP"], marker="o", markersize=3, linestyle=linestyle,  uplims=True, lolims=True, lw=lw, elinewidth = 0.5, label= "DFT-pbe Reference", c="orange", capsize=0.5)
   
    
    if reference is not None:
        axs[0].plot(reference[:,0], reference[:,1], "--o", markersize = 3, linewidth=lw, label = "Literature", c="green")
        
    axs[0].set_xlabel("Density $[g/cm^3]$")
    axs[0].set_ylabel("Pressure [GPa]")
    axs[0].set_title(f"Pressure at Different Densities ")
    axs[0].legend()

    if zoomsize is not None:
        axins = axs[0].inset_axes(zoomsize)
        axins.errorbar(x=x, y = df["LP"], yerr=df["sigma_LP"], marker="o", markersize=3, linestyle = linestyle, lw = lw,  uplims=True, lolims=True,  elinewidth=0.5, label ="GAP",c="blue",capsize=0.5)
        axins.errorbar(x=xRef, y = dfRef["RP"], yerr=dfRef["sigma_RP"], marker="o", markersize=3, linestyle = linestyle, lw=lw, uplims=True, lolims=True, elinewidth = 0.5, label= "DFT-pbe Reference", c="orange",capsize=0.5)
        #if df_noRef is not None:
        #    axins.errorbar(x=xnoRef, y = df_noRef["LP"], yerr=df_noRef["sigma_LP"], marker="o",  markersize=3, linestyle=linestyle, lw = lw, uplims=True, lolims=True, elinewidth=0.5,c="blue", capsize=0.5)
        if xlimit is not None:
            axins.set_xlim(xlimit)
        if ylimit is not None:
            axins.set_ylim(ylimit)
        if reference is not None:
            axins.plot(reference[:,0], reference[:,1], "--o",  markersize = 3, linewidth=lw, c="green")

            
    axs[1].errorbar(x=x, y = df["LE"], yerr=df["sigma_LE"], marker="o", markersize=3, linestyle=linestyle, uplims=True, lolims=True, lw=lw, elinewidth=0.5, label ="GAP", c="blue", capsize=0.5)
    axs[1].errorbar(x=xRef, y = dfRef["RE"], yerr=dfRef["sigma_RE"], marker="o", markersize=3,  linestyle=linestyle,uplims=True, lolims=True, lw=lw, elinewidth=0.5, label= "DFT-pbe Reference",  c="orange", capsize=0.5)
    #if df_noRef is not None:
    #    axs[1].errorbar(x=xnoRef, y = df_noRef["LE"], yerr=df_noRef["sigma_LE"], marker="o", markersize=3, uplims=True, lolims=True, lw=0, elinewidth=0.5, c="blue", capsize=0.5)

    axs[1].set_xlabel("Density $[g/cm^3]$")
    axs[1].set_ylabel("Energy/atom [Hartrees]")
    axs[1].set_title("Energy at Different Densities")
    axs[1].legend()
    
    plt.suptitle(f"EOS at {T} K")
