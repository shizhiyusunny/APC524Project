#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from . import Postprocessing
from fenics import *

class toNPY(Postprocessing.Postprocessing):
    def store(t,f):
        u=f.displacement
        arr_u=u.compute_vertex_values()
        fileD = File("data/displacement.txt")
        fileD << arr_u

