#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from abc import ABC, abstractmethod

class Postprocessing(ABC):
    @abstractmethod
    def store(t,f):
        pass
    def postprocess(t,process,flag):
        from . import toPVD,toNPY,toPLT
        for i in process:
            if i=="pvd" and flag==0:
                P=toPVD.toPVD()
                P.store(t,f)
            elif i=="npy" and flag==0:
                P=toNPY.toNPY()
                P.store(t,f)
            elif i=="pltDisp" and flag==1:
                P=toPLT.toDisPLT()
                P.store(t,f)

