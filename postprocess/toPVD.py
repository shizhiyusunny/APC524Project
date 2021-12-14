#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from . import Postprocessing
from fenics import *

class toPVD(Postprocessing.Postprocessing):
    def store(t,f):
        u=f.displacement
        u.rename("displacement","")
        fileD = File("data/displacement.pvd")
        fileD << u

