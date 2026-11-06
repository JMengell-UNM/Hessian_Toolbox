'''
Migrating Runtime Test Code to Separate Script
'''

#%%
#
from Hessian_Toolbox import Hessian

H2O = Hessian('./Test_Hessians/H2O_vib.hess')


#%%
#Testing Normal Mode Comination
from Hessian_Toolbox import Hessian
benzene = Hessian('./Test_Hessians/benzene_Zindo_Vib.hess')

benzene.gen_phase_displaced_oscilations(6,7,8)
