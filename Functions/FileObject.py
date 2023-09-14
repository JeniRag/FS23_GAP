# -*- coding: utf-8 -*-
"""
Created on Sun Mar 26 19:43:45 2023

@author: Sujeni
"""
import numpy as np
import re

class FileObject():
    """
    A class to read the output files of a NVT simulation. Customized to my specific applications.
    """
    def __init__(self, file, file_type, dt=None, Temperature=None, N = None):
        """
        

        Parameters
        ----------
        file : str
            Path to File.
        file_type : str
            Type of File. Available options: "QE" or "i-Pi.out"
        dt : float 
            Timestep of Simulation in a.u. The default is None.
        Temperature : int
            Temperature of Simulation. The default is None.
        N : int 
            Number of Atoms in the simulation. The default is None.

        Returns
        -------
        None.

        """
        
        self.file =file
        self.file_type = file_type
        
        self.time = None
        self.steps=None
        self.Energy=None
        self.Pressure_tot=None
        self.Pressure_elec=None
        self.Volume=None
        self.Temperature=Temperature
        self.KE=None
        self.PE=None
        self.dt =dt
        self.label=None
        self.N=N #number of atoms
        
        
        self.extractInfo()
        
        
        
    def extractInfo(self):

        #conversion factors
        bar_to_GPa = 1e-4
        au_to_s = 4.83781*1e-17 #2.4189* 1e-17
        s_to_ps = 1e12 # seconds to pico second
        Ry_to_ev =  13.605
        kB = 8.617*1e-5 #Boltzmann constant in eV/K
        
        
        if self.file_type=="QE":   
   
            if self.dt==None:
                raise ValueError(
                    f"dt not provided."
                )
                
            if self.Temperature==None:
                raise ValueError(
                    f"Temperature not provided."
                )
                
            
            
            f=open(self.file, "r")
            lines=f.readlines()[1:]
    
            Nframes = len(lines)
            
            #empty arrays to store values
            steps = np.zeros((Nframes, 1), dtype=float) #iteration
            Energy = np.zeros((Nframes, 1), dtype=float)
            Pressure_elec = np.zeros((Nframes, 1), dtype=float)
            Pressure_tot = np.zeros((Nframes, 1), dtype=float)
        
            #extract and store values    
            for i,x in enumerate(lines): 
                Energy[i] = float(x.split()[0]) #Ry Energy per atom
                Pressure_elec[i] = float(x.split()[1])
                Pressure_tot[i] = float(x.split()[2])
                steps[i] = i #float(x.split()[3]) #columns
        
            self.Pressure_tot=Pressure_tot #GPa
            self.Energy_elec=Energy * Ry_to_ev  # in eV
            self.Energy = self.Energy_elec + 1/2 * kB * self.Temperature  #self.N/2 * kB * self.Temperature #total Energy
            self.Pressure_elec = Pressure_elec # in GPa
            self.steps=steps
            self.time = self.dt *steps* au_to_s * s_to_ps
            self.dt = self.dt * au_to_s * s_to_ps # convert a.u to ps
        
                
            f.close()
            
        elif self.file_type=="i-Pi.out":
        # i-Pi trajectory output
        
            f=open(self.file, "r")
            lines = f.readlines()
            headers = lines[:12]
            values = lines[12:]
            
            #detect current units of file
            for h in headers:
                match_pressure = re.search(r'pressure_md\{(.+?)\}', h)
                if match_pressure:
                    pressure_units = match_pressure.group(1)
                
                match_time = re.search(r'time\{(.+?)\}', h)
                if match_time:
                    time_units = match_time.group(1)
                    
                match_KE= re.search(r'kinetic_md\{(.+?)\}', h)
                if match_KE:
                    KE_units = match_KE.group(1)
                
                match_PE= re.search(r'potential\{(.+?)\}', h)
                if match_PE:
                    PE_units = match_PE.group(1)
                    
                
                    
           
            Nframes = len(values)
            
            #empty arrays to store values
            steps = np.zeros((Nframes, 1), dtype=float)
            time = np.zeros((Nframes, 1), dtype=float)
            conservedE = np.zeros((Nframes, 1), dtype=float)
            KE = np.zeros((Nframes, 1), dtype=float)
            PE = np.zeros((Nframes, 1), dtype=float)
            Pressure = np.zeros((Nframes, 1), dtype=float)
            Volume = np.zeros((Nframes, 1), dtype=float)
            Temperature = np.zeros((Nframes, 1), dtype=float)
        
            i = 0
            for x in values: 
                steps[i] = float(x.split()[0]) #columns
                time[i] = float(x.split()[1]) #pico second
                conservedE[i]= float(x.split()[2]) #eV
                Temperature[i] = float(x.split()[3])
                KE[i] = float(x.split()[4]) # The kinetic energy of the (extended) classical system. average in eV
                PE[i] = float(x.split()[5]) # The pressure of the (extended) classical system.
                Pressure[i] = float(x.split()[6]) #bar
                Volume[i] = float(x.split()[7])
                i += 1
                
            self.steps = steps
            self.time = time
            self.conservedE = conservedE
            self.Temperature = Temperature
            self.KE = KE
            self.PE = PE
            self.Pressure_tot = Pressure
            
            #convert to desired unit
            if pressure_units=="bar":
                self.Pressure_tot *= bar_to_GPa
                
            elif pressure_units=="GPa" or pressure_units == "gigapascal":
                pass
       
            else:
                raise ValueError(
                    f"Pressure unit {pressure_units} not implemented or detected.")
                
                
            if time_units=="picosecond":
                pass
            else:
                raise ValueError(
                    f"Time unit {time_units} not implemented or detected.")
            
            if KE_units=="electronvolt":
                pass
            else:
                raise ValueError(
                   f"Energy unit {KE_units} not implemented or detected.")
              
            if PE_units=="electronvolt":
                 pass
            else:
                 raise ValueError(
                    f"Energy unit {KE_units} not implemented or detected.")
                
            
            self.Volume = Volume
            self.Energy = KE+PE #eV
            self.dt =self.time[1] - self.time[0]
            
           
            f.close()
         
        else:
            raise ValueError(
                f"File Type not supported. Possible Formats "
                "are 'i-Pi.out', 'QE'. "
            )
            
        return
    
    