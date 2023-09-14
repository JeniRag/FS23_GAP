# -*- coding: utf-8 -*-
"""
Created on Sun Mar 26 19:43:45 2023

@author: Sujeni
"""

class FileObject():
    def __init__(file, file_type):
        self.file =file
        self.file_type = file_type
        self.time = None
        self.steps=None
        self.Energy=None
        self.Pressure_tot=None
        self.Pressure_elec=None
        self.Volume=None
        self.Temperature=None
        self.KE=None
        self.PE=None
        
        self.extractInfo()
        
        
        
    def extractInfo(self):
        """
        Parameters
        ----------
        file : file path
        file_type : LAMMPS or QE
    
        Returns
        -------
        steps, time, Energy, Pressure
    
        """
        #conversion factors
        bar_to_GPa = 1e-4
        au_to_s = 2.4189* 1e-17
        s_to_ps = 1e12 # seconds to pico second
        
        
        if self.file_type=="QE":
            f=open(self.file, "r")
            lines=f.readlines()[1:]
    
            Nframes = len(lines)
            steps = np.zeros((Nframes, 1), dtype=float) #iteration
            Energy = np.zeros((Nframes, 1), dtype=float)
            Pressure_elec = np.zeros((Nframes, 1), dtype=float)
            Pressure_tot = np.zeros((Nframes, 1), dtype=float)
        
            i=0
            for x in lines: 
                Energy[i] = float(x.split()[0])
                Pressure_elec[i] = float(x.split()[1])
                Pressure_tot[i] = float(x.split()[2])
                steps[i] = float(x.split()[3]) #columns
        
                i+=1
            
            self.Pressure_tot=Pressure
            self.Energy=Energy
            self.Pressure_elec = Pressure_elec
            self.steps=steps
        
                
            f.close()
            
        elif file_type=="LAMMPS":
            f=open(file, "r")
            lines=f.readlines()[12:]
        
            Nframes = len(lines)
            steps = np.zeros((Nframes, 1), dtype=float)
            time = np.zeros((Nframes, 1), dtype=float)
            KE = np.zeros((Nframes, 1), dtype=float)
            PE = np.zeros((Nframes, 1), dtype=float)
            Pressure = np.zeros((Nframes, 1), dtype=float)
            Volume = np.zeros((Nframes, 1), dtype=float)
            Temperature = np.zeros((Nframes, 1), dtype=float)
        
            i = 0
            for x in lines: 
                steps[i] = float(x.split()[0]) #columns
                time[i] = float(x.split()[1]) #pico second
                Temperature[i] = float(x.split()[3])
                KE[i] = float(x.split()[4]) # The kinetic energy of the (extended) classical system.
                PE[i] = float(x.split()[5]) # The pressure of the (extended) classical system.
                Pressure[i] = float(x.split()[6]) #bar
                Volume[i] = float(x.split()[7])
                i += 1
                
            self.steps=steps
            self.time = time
            self.Temperature =Temperature
            self.KE = KE
            self.PE=PE
            self.Pressure_tot =Pressure * bar_to_GPa
            self.Volume=Volume
            self.Temperature
            self.Energy = KE #look if KE or PE??
            
           
            f.close()
         
        else:
            print("file type not implemented")
            
        return
    
    