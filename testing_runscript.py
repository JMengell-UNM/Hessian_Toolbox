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

benzene.combine_normal_modes(6,7,8)
