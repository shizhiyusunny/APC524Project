from fenics import *
class FunctionW:
    def __init__(self, inputs, mesh, lm):
        degreeElements = 1
        self.mesh = mesh
        self.dim = inputs['dim']
        if lm == []:
            MFS = VectorFunctionSpace(mesh, 1, inputs['dim'])
            u  = Function(MFS)
            self.displacement = u
            self.unknown = u
            return
        P1 = FiniteElement('Lagrange', mesh.ufl_cell(), degreeElements)
        if inputs['dim'] == 1:
            disp_space = P1
        elif inputs['dim'] == 2:
            disp_space = (P1*P1)
        elif inputs['dim'] == 3:
            disp_space = (P1*P1*P1)
        mixed_space = [disp_space]
        R = FiniteElement('Real', mesh.ufl_cell(), 0)
        lm_dofs = []
        for l in lm:
            if l.dof == 'translation':
                if inputs['dim'] == 1:
                    tr_space = R
                elif inputs['dim'] == 2:
                    tr_space = (R*R)
                elif inputs['dim'] == 3:
                    tr_space = (R*R*R)
                mixed_space.append(tr_space)
                lm_dofs.append('translation')
            if l.dof == 'rotation':
                if inputs['dim'] == 2:
                    rot_space = R
                elif inputs['dim'] == 3:
                    rot_space = (R*R*R)
                mixed_space.append(rot_space)
                lm_dofs.append('rotation')
        MFS = FunctionSpace(mesh, MixedElement(mixed_space))
        f  = Function(MFS)
        if 'translation' not in lm_dofs:
            u, c_rot = split(f)
            self.lm_rot = c_rot
        elif 'rotation' not in  lm_dofs:
            u, c_trans = split(f)
            self.lm_trans = c_trans
        else:
            u, c_trans, c_rot = split(f)
            self.lm_trans = c_trans
            self.lm_rot = c_rot
        self.displacement = u
        self.unknown = f
