import numpy as np
import yaml
class Optimizer:
    def builder(inputs):
        if 'Optimizer' in inputs:
            op = Optimizer(inputs)
            return op
        else:
            return EmptyOptimizer()

    def __init__(self, inputs):
        inp = inputs['Optimizer']
        self.divisions = np.geomspace(inp['Min Divisions'],inp['Max Divisions'],num=10)
        self.tolerance = inp['Tolerance']
        self.first = True
        inputs['Mesh']['Division'] = self.divisions[0]
        print('Divisions set to: ', inputs['Mesh']['Division'])
        self.inputs = inputs
        self.div_index = 0
        self.not_converged = True

    def check(self, func):
        new_val = func.calc_magnitude()
        if self.first:
            self.first = False
        else:
            if abs(self.val - new_val) < self.tolerance:
                self.not_converged = False
                return self.not_converged
        self.val = new_val
        self.div_index += 1
        if self.div_index == np.shape(self.divisions)[0] - 1:
            return False
        self.inputs['Mesh']['Division'] = self.divisions[self.div_index]
        print('Divisions set to: ', self.inputs['Mesh']['Division'])
        return self.not_converged

class EmptyOptimizer:
    def __init__(self):
        pass
    def check(self, func):
        return False
