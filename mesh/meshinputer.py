from fenics import *
from mshr import *
from dolfin import *
import filetype
import os


class MeshInputer:

    def __init__(self, filename):
        self.filename = filename
        curpath = os.path.dirname(os.path.realpath("test.ipynb"))
        self.yamlpath = os.path.join(curpath, self.filename)

    def get_extension(self):
        kind = filetype.guess(self.filename)
        if kind is None:
            print("Cannot guess file type!")
            return
        return kind.extension

    def mesh_input(self):
        ext = self.get_extension()
        ext_list = ['ele', '.node', '.mesh', '.msh', 'gmsh', 'grid', 'inp', 'e', 'exo', 'ncdf', 'vrt', 'cell']

        if ext == 'xml':
            mesh = Mesh(self.yamlpath)
        elif ext in ext_list:
            print("This meshfile is not in .xml format. Use dolfin-convert in bash to get .xml file.")
            return
        else:
            print("Unfortunately, the format of this meshfile does not our requirements.")
            return

        return mesh
