#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
from . import Postprocessing
from . import ymlReader

class Annealer():
    def ConstantStep(res,f,inputs):
        n, process=ymlReader.Reader.read(inputs)
        arr_step=np.linspace(0,1,n)
        flag_prog=0
        for t in arr_step:
            #res.updatetraction(t)
            res.calculate_residual(f)
            res.solve()
            Postprocessing.Postprocessing.postprocess(t,process,flag_prog)
        flag_prog=1
        Postprocessing.Postprocessing.postprocess(t,process,flag_prog)
            
            
    

