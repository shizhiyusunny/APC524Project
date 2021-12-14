#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from fenics import *
class Reader():
    def read(inputs):
        Nsteps=0
        process=[]
        for key, val in inputs.items():
            if key == 'AnnealingStep':
                Nsteps=val
            if key == 'process':
                for v in val:
                    process.append(v)
        return Nsteps, process

