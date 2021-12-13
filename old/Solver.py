#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from fenics import *
import numpy as np
import constants as cst
class Solver():
    def __init__(self,mesh):
        self.mesh=mesh
        self.tf=cst.tfinal
        self.delt=cst.tstep
        u=[]   
        c_trans=[]
        c_rot=[]
        strain=[]
        stress=[]
        disp=[]
    def constantStep(self):
        arr_time=np.linspace(0,self.tf,self.delt)
        ###Append Initial Conditions Here!!!
        for i in arr_time:
            Res=Energyfunction() ####
            mySim=Simulator(self.mesh)
            u_i,c_trans_i,c_rot_i, strain_i, stress_i, disp_i=mySim.stepForward(Res)
            self.u.append(u_i)
            self.c_trans.append(c_trans_i)
            self.c_rot.append(c_rot_i)
            self.strain.append(strain_i)
            self.stress.append(stress_i)
            self.disp.append(disp_i)
        return self.u, self.strain, self.stress, self.disp
    def adaptiveStep(self):
        pass

        

