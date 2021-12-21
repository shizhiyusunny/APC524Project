from . import annealer
import numpy as np

class ConstantStep(annealer.Annealer):
    def __init__(self, inputs):
        self.n = inputs['Number of Steps']
    def stepper(self, res, f, processes = []):
        arr_step = np.linspace(0,1,self.n + 1)
        for t in arr_step:
            res.update_tractions(float(t))
            res.calculate_residual(f)
            res.solve()
            annealer.Annealer.store_all(processes, f, t)
