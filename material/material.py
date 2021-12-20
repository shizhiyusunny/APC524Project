from fenics import *
from dolfin import *

class MaterialConstant(UserExpression):
    def builder(markers, inputs):
        from . import heterogeneous
        material_constants = {}
        for key, val in inputs.items():
            if isinstance(val, list):
                material_constants[key] = heterogeneous.HeterogeneousConstant(markers, val)
            else:
                material_constants[key] = val
        MaterialConstant.processor(material_constants)
        return material_constants
    def processor(material_constants):
        if 'E' in material_constants and 'nu' in material_constants:
            E = material_constants['E']
            nu = material_constants['nu']
            material_constants['mu'] = E/2/(1+nu)
            material_constants['lambda'] = E*nu/(1-nu*nu)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
    def eval_cell(self, value, x, ufc_cell):
        pass
    def value_shape(self):
        return ()
