# %%

import numpy as np
import matplotlib.pyplot as plt
import os

class Hessian:
    '''Adapted from Indu's LLM generated code
    '''
    def __init__(self, file_path):
        self.file_path=file_path
        self.hess_dict={}
        
        with open(file_path, 'r') as file:
            contents=file.readlines()
            i=0
            while i<len(contents):
                
                if contents[i].strip()=='$vibrational_frequencies':
                    j=2
                    freqs=[]
                    while contents[i+j].strip()!='':
                        freqs.append(float(contents[i+j].split()[1]))
                        j+=1
                        self.hess_dict['vibrational_frequencies']=freqs
                    del freqs
                
                if contents[i].strip()=='$normal_modes':
                    i+=1 #next line
                    num_rows, num_cols =[int(n) for n in contents[i].strip().split()]
                    i+=1
                    normal_modes=np.zeros((num_rows, num_cols))
                    col=0
                    while col<num_cols:
                        i+=1
                        block_size=min(5, num_cols-col)
                        for row in range(num_rows):
                            line = contents[i].split()
                            for k in range(block_size):
                                normal_modes[row, col + k] = float(line[1 + k])
                            i += 1
                        col += block_size
                                        
                    self.hess_dict['normal_modes']=normal_modes
                    del num_cols,normal_modes
                    
                if contents[i].strip()=='$atoms':
                    self.hess_dict['atomic_symbols']=[]
                    self.hess_dict['molecular_weights']=[]
                    #Formating atomic coords to match Q
                    self.hess_dict['atomic_coords']=np.zeros(num_rows)
                    
                    
                    j=2
                    k=0
                    while contents[i+j].strip()!='':
                        symbol, MW, x, y, z = contents[i+j].split()
                        
                        self.hess_dict['atomic_symbols'].append(symbol)
                        self.hess_dict['molecular_weights'].append(float(MW))
                        for n in [x, y, z]:
                            self.hess_dict['atomic_coords'][k]=n
                            k+=1
                                                                
                        j+=1
                        
                    del symbol, MW, x, y, z, num_rows
                    
                if contents[i].strip()=='$frequency_scale_factor':
                    self.hess_dict['frequency_scale_factor']=float(contents[i+1].strip())
              
                i+=1
            del contents


    def step_xyz(self, Qn, step_size):
        ''' 
        Qn is the index value of the vibrational mode being called. 
        Step size is the density of steps between zero and max distortion. Value between 0 and 1.
        '''
        Mw = self.hess_dict['molecular_weights']
        R = self.hess_dict['atomic_coords']
        Q = self.hess_dict['normal_modes'][:,Qn]

        Mw_xyz = np.repeat(np.array(Mw), 3)
        #print(*[Mw_xyz, Q, R, step_size],sep='\n\n')
        Q_displaced = step_size * Q / np.sqrt(Mw_xyz) + R # Calculates a displaced geometry, Mw weighted
        # Q_displaced = Q_displaced.reshape(-1,3) # rearranges matrix into rows of [x y z] for input into .xyz
        return Q_displaced

    def write_xyz(self, filename, geometry, comment="Generated via VSOC Python Script"):
        n_atoms = len(self.hess_dict['atomic_symbols'])
        geometry = (np.asarray(geometry)).reshape(-1,3) # Rearranges matrix in rows of [x y z] for input in to .xyz
        symbols = self.hess_dict['atomic_symbols']

        script_dir = os.path.dirname(os.path.abspath(__file__))
        filepath = os.path.join(script_dir, filename + ".xyz") # Defines the location of the .py script as where the .xyz files should be generated. Change to working_directory for alternative behavior

        with open(filepath, "w") as f:
            f.write(f"{n_atoms}\n")
            f.write(f"{comment}\n")
            for symbol, (x, y, z) in zip(symbols, geometry):
                 f.write(f"{symbol:<2s}  {x:>15.10f}  {y:>15.10f}  {z:>15.10f}\n")

    def write_inp(self, filename, charge=0, mult=1, xc_functional="B3LYP", basis_set="ZORA-def2-SVP"):
        script_dir = os.path.dirname(os.path.abspath(__file__))
        filepath = os.path.join(script_dir, filename + ".inp")

        with open(filepath, "w") as f:
            f.write(f"! {xc_functional} {basis_set} SARC/J VerySlowConv noTRAH LargePrint PrintBasis\n\n")

            f.write("%basis\n")
            f.write('\tnewGTO Pt "SARC-ZORA-TZVP" end\n') #maybe parameterize the seperate Pt basis
            f.write("\tend\n\n")

            f.write("%scf\n")
            f.write("\tmaxiter 300\n")
            f.write("\tend\n\n")
        
            f.write("%tddft\n")
            f.write("\tnroots 3\n")
            f.write("\ttriplets True\n")
            f.write("\tdoSOC True\n")
            f.write("\tend\n\n")

            f.write("%rel\n")
            f.write("\t#Default values\n")
            f.write("\tSOCType 3\n")
            f.write("\tSOCFlags 1,4,3,0\n")
            f.write("\tend\n\n")

            f.write(f"*xyzfile {charge} {mult} {filename}.xyz\n")

    def gen_normal_mode_xyz(self, *Qns, n_steps=None, max_displacement=1):
        step_size = max_displacement / n_steps
        step_list = np.linspace(0, max_displacement, num=n_steps, endpoint=True)
        for Qi in Qns:
            # Q = self.hess_dict['normal_modes'][:,Qi] # unnecesssary geoemtry definition
            for step in step_list:
                Title = f"{step}_step_displacement_vib_mode_{Qi}" # Iteriative file title for the given mode and step
                Q_displaced = self.step_xyz(Qi, step_size*step) # Generates each geometry for each vib and each step
                self.write_xyz(filename=Title, geometry=Q_displaced)
                #self.write_inp(filename=Title)
                ''' iterate this such that all .inp files are written simultainiously'''
                
        '''
        for Qn in Qns:
            displacement_grid= [ 'a function of stepsize and max displacement']

            for D in displacement_grid:
                displaced_geom = self.step_xyz(Qn, D)
                self.write_xyz(displaced_geom, target_dir)
        
        '''

    def gen_normal_mode_inputs(self, *Qns, n_steps=None, max_displacement=1):
        step_list = np.linspace(0, max_displacement, num=n_steps, endpoint=True)
        for Qi in Qns:
            for step in step_list:
                Title = f"{step}_step_displacement_vib_mode_{Qi}" # Iteriative file title for the given mode and step
                self.write_inp(filename=Title)

    def gen_xyz_inp(self, *Qns, n_steps=None, max_displacement=1):
        self.gen_normal_mode_xyz(*Qns, n_steps=n_steps, max_displacement=max_displacement)
        self.gen_normal_mode_inputs(*Qns, n_steps=n_steps, max_displacement=max_displacement)

        '''WORKS'''

    def combine_freq(self, *Qns, ):
         ''' 
        Qn is the index value of the vibrational mode being called. 
        '''
        # Mw = self.hess_dict['molecular_weights']
        # R = self.hess_dict['atomic_coords']
        # for Qn in Qns:
        #     Q = self.hess_dict['normal_modes'][:,Qn]
        
        # Mw_xyz = np.repeat(np.array(Mw), 3)
    pass

    def write_hess():
        pass

    def dos_histplot(self, bin_size=50):
        all_freq = self.hess_dict['vibrational_frequencies']
        all_freq = np.array(all_freq, dtype=float).flatten()

        min_value = int(np.floor(all_freq.min()))
        max_value = int(np.ceil(all_freq.max()))
        custom_bins = np.arange(min_value, max_value + bin_size, bin_size)

        plt.hist(all_freq, bins=custom_bins)

        plt.title(f"Density of Vibrational States, Bin = {bin_size}", fontsize=14)
        plt.xlabel(r'Wavenumber ($cm^{-1}$)', fontsize=12)
        plt.ylabel("Number of States", fontsize=12)
        plt.show()

# %%
''' All Code below is for testing purposes'''

H2 = Hessian(r"C:\Users\apelr\OneDrive - University of New Mexico\Python Scripts\Josh_Scripts\H2_vib.hess")
# BDT = Hessian(r"C:\Users\apelr\OneDrive - University of New Mexico\Python Scripts\Josh_Scripts\BpyPtBDT_fine_vib.hess")

script_dir = os.path.dirname(os.path.abspath(__file__))
# Combine the directory with your filename
file_path = os.path.join(script_dir, "BpyPtBDT_fine_vib.hess")
# Pass the dynamic path to your Hessian function
BDT = Hessian(file_path)

                
# Calculate the displaced geometry

# script testing
# print(H2.step_xyz(-1, 1))
# print(H2.hess_dict['normal_modes'])
# H2.write_xyz(filename="test_H2", geometry=H2.step_xyz(-1, 1), )

# print(BDT.step_xyz(-1, 1))
# print(BDT.hess_dict['normal_modes'])
# print(BDT.hess_dict['atomic_coords'].reshape(-1,3))
# print(BDT.hess_dict['atomic_symbols'])
# BDT.write_xyz(filename="test_BDT", geometry=BDT.step_xyz(14, 2))
# print(BDT.hess_dict['vibrational_frequencies'])

# BDT.gen_xyz_inp(5, 7, n_steps=3, max_displacement=3)

BDT.dos_histplot(bin_size=20)

            
