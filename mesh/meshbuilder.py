from fenics import *
from mshr import *
from dolfin import *


class MeshBuilder:

    def __init__(self, divisions):
        self.divisions = divisions

    # Use generate_mesh to build mesh, first you need a domain.
    # The cell size is about the characteristic length over self.divisions.
    def Generate(self, domain):
        mesh = generate_mesh(domain, *self.divisions)
        return mesh

    # Encapsulate different dimensions,
    # return the mesh for a d-dimensional cube given a list or tuple
    # specifying the division into cells in the spatial coordinates.
    def UnitHyperCube(self):
        mesh_classes = [UnitIntervalMesh, UnitSquareMesh, UnitCubeMesh]
        d = len(self.divisions)
        mesh = mesh_classes[d - 1](*self.divisions)
        return mesh

    def HyperCube(self, length):
        mesh = self.UnitHyperCube()
        mesh.coordinates()[:] = length * mesh.coordinates()[:]
        return mesh

    def MeshRect(self, center, length, breadth):
        centerx = center[0]
        centery = center[1]
        pt1x = centerx - 0.5 * length
        pt1y = centery - 0.5 * breadth
        pt2x = centerx + 0.5 * length
        pt2y = centery + 0.5 * breadth
        mesh = RectangleMesh(Point(pt1x, pt1y), Point(pt2x, pt2y), *self.divisions)
        return mesh

    def MeshBox(self, center, length, width, height):
        centerx = center[0]
        centery = center[1]
        centerz = center[2]
        pt1x = centerx - 0.5 * length
        pt1y = centery - 0.5 * width
        pt1z = centerz - 0.5 * height
        pt2x = centerx + 0.5 * length
        pt2y = centery + 0.5 * width
        pt2z = centerz + 0.5 * height
        mesh = BoxMesh(Point(pt1x, pt1y, pt1z), Point(pt2x, pt2y, pt2z), *self.divisions)
        return mesh

    def MeshUnitCircle(self):
        mesh = UnitCircleMesh(*self.divisions)
        return mesh

    def MeshUnitSphere(self):
        mesh = UnitSphereMesh(*self.divisions)
        return mesh

    def MeshCircle(self, r):
        mesh = self.MeshUnitCircle()
        mesh.coordinates()[:] = r * mesh.coordinates()[:]
        return mesh

    def MeshSphere(self, r):
        mesh = self.MeshUnitSphere()
        mesh.coordinates()[:] = r * mesh.coordinates()[:]
        return mesh

