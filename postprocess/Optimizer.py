#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from fenics import *
import numpy as np
import constants as cst
class Optimizer():
    def __init__(self,mesh):
        self.eps=cst.threshold #read from constants from user input
        self.maxmesh=cst.meshsize #read from constants from user input
        self.iterations=cst.iterations
        self.error=0
        self.optsize=0
        self.mesh=mesh

    def optimizeMesh(self):
        while self.error>=self.eps:
            if self.optsize>self.maxmesh:
                break
            else:
                arr_mesh=np.linspace(0,self.maxmesh,self.iterations).astype(int)
                disp_last=0
                for i in arr_mesh:
                    #Remash the whole structure
                    #Solve it
                    mySolver=Solver(myMesh)
                    res=mySolver.constantStep()
                    disp_now=np.average(res[3][-1][-1])
                    error=(disp_now-disp_last)/disp_last
                    disp_last=dis_now
                    self.optsize=i
        return self.optsize
        
                
                    
                    
                

