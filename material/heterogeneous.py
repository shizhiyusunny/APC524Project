from . import material
from fenics import *
from dolfin import *

class HeterogeneousConstant(UserExpression):
    def __init__(self, markers, val, **kwargs):
        self.val = {}
        for v in val:
            self.val[v['marker']] = v['value']
        self.markers = markers
        super().__init__(**kwargs)

    def eval_cell(self, value, x, ufc_cell):
        value[0] = self.val[self.markers[ufc_cell.index]]

    def value_shape(self):
        return ()
