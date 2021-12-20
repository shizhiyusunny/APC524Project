from fenics import *
from mshr import *
from dolfin import *
from .shape import Shape
from .setshape import SetShape
from .meshinputer import MeshInputer
from .meshbuilder import MeshBuilder


class MeshGenerator:

    def __init__(self, meshdict):
        self.meshdict = meshdict

    division_1d = [20]
    division_2d = [20, 20]
    division_3d = [20, 20, 20]

    @staticmethod
    def get_list1d(n):
        list_n = [n]
        return list_n

    @staticmethod
    def get_list2d(str_n):
        list_n = []
        split_n = str_n.split(',')
        list_nx = float(split_n[0])
        list_ny = float(split_n[1].lstrip())
        list_n.append(list_nx)
        list_n.append(list_ny)
        return list_n

    @staticmethod
    def get_list3d(str_n):
        list_n = []
        split_n = str_n.split(',')
        list_nx = float(split_n[0])
        list_ny = float(split_n[1].lstrip())
        list_nz = float(split_n[2].lstrip())
        list_n.append(list_nx)
        list_n.append(list_ny)
        list_n.append(list_nz)
        return list_n

    def get_shape2d(self, shape_dict):
        str_center = shape_dict['Center']
        center = self.get_list2d(str_center)
        shape = Shape(center)
        return shape

    def get_shape3d(self, shape_dict):
        str_center = shape_dict['Center']
        center = self.get_list3d(str_center)
        shape = Shape(center)
        return shape

    def domain_generator(self, key, shape_dict):
        domain = None

        if key == 'Interval':
            center = shape_dict['Center']
            shape = Shape(center)
            domain = shape.interval(shape_dict['Length'])

        elif key == 'Circle':
            shape = self.get_shape2d(shape_dict)
            domain = shape.circle(shape_dict['Radius'])

        elif key == 'Rectangle':
            shape = self.get_shape2d(shape_dict)
            domain = shape.rectangle(shape_dict['Length'], shape_dict['Breadth'])

        elif key == 'Square':
            shape = self.get_shape2d(shape_dict)
            domain = shape.square(shape_dict['Length'])

        elif key == 'Sphere':
            shape = self.get_shape3d(shape_dict)
            domain = shape.sphere(shape_dict['Radius'])

        elif key == 'Box':
            shape = self.get_shape3d(shape_dict)
            domain = shape.box(shape_dict['Length'], shape_dict['Width'], shape_dict['Height'])

        elif key == 'Cube':
            shape = self.get_shape3d(shape_dict)
            domain = shape.cube(shape_dict['Length'])

        return domain

    def domain_operation(self):
        key = list(self.meshdict.items())[0][0]
        domain1 = self.domain_generator(key, self.meshdict[key])

        for key, val in self.meshdict.items():
            if key == 'Subtract':
                for key_sub, val_sub in self.meshdict['Subtract'].items():
                    domain2 = self.domain_generator(key_sub, self.meshdict['Subtract'][key_sub])
                    setshape = SetShape(domain1, domain2)
                    domain1 = setshape.subtraction()

            elif key == 'Add':
                for key_add, val_add in self.meshdict['Add'].items():
                    domain2 = self.domain_generator(key_add, self.meshdict['Add'][key_add])
                    setshape = SetShape(domain1, domain2)
                    domain1 = setshape.addiction()

        return domain1

    def mesh_builder(self):
        if 'Method' in self.meshdict:
            if self.meshdict['Method'] == 'generate_mesh':
                domain = self.domain_operation()
                if 'Division' in self.meshdict:
                    division = self.get_list1d(self.meshdict['Division'])
                else:
                    division = self.division_1d
                meshbuilder = MeshBuilder(division)
                mesh = meshbuilder.Generate(domain)

            if self.meshdict['Method'] == 'UnitIntervalMesh':
                if 'Division' in self.meshdict:
                    division = self.get_list1d(self.meshdict['Division'])
                else:
                    division = self.division_1d
                meshbuilder = MeshBuilder(division)
                mesh = meshbuilder.UnitHyperCube()

            if self.meshdict['Method'] == 'UnitSquareMesh':
                if 'Division' in self.meshdict:
                    division = self.get_list2d(self.meshdict['Division'])
                else:
                    division = self.division_2d
                meshbuilder = MeshBuilder(division)
                mesh = meshbuilder.UnitHyperCube()

            if self.meshdict['Method'] == 'UnitCubeMesh':
                if 'Division' in self.meshdict:
                    division = self.get_list3d(self.meshdict['Division'])
                else:
                    division = self.division_3d
                meshbuilder = MeshBuilder(division)
                mesh = meshbuilder.UnitHyperCube()

            if self.meshdict['Method'] == 'IntervalMesh':
                if 'Division' in self.meshdict:
                    division = self.get_list1d(self.meshdict['Division'])
                else:
                    division = self.division_1d
                meshbuilder = MeshBuilder(division)
                mesh = meshbuilder.HyperCube(self.meshdict['Interval']['Length'])

            if self.meshdict['Method'] == 'SquareMesh':
                if 'Division' in self.meshdict:
                    division = self.get_list2d(self.meshdict['Division'])
                else:
                    division = self.division_2d
                meshbuilder = MeshBuilder(division)
                mesh = meshbuilder.HyperCube(self.meshdict['Square']['Length'])

            if self.meshdict['Method'] == 'CubeMesh':
                if 'Division' in self.meshdict:
                    division = self.get_list3d(self.meshdict['Division'])
                else:
                    division = self.division_3d
                meshbuilder = MeshBuilder(division)
                mesh = meshbuilder.HyperCube(self.meshdict['Cube']['Length'])

            if self.meshdict['Method'] == 'RectangleMesh':
                if 'Division' in self.meshdict:
                    division = self.get_list2d(self.meshdict['Division'])
                else:
                    division = self.division_2d
                meshbuilder = MeshBuilder(division)
                center = self.get_list2d(self.meshdict['Rectangle']['Center'])
                mesh = meshbuilder.MeshRect(center, self.meshdict['Rectangle']['Length'],
                                            self.meshdict['Rectangle']['Breadth'])

            if self.meshdict['Method'] == 'BoxMesh':
                if 'Division' in self.meshdict:
                    division = self.get_list3d(self.meshdict['Division'])
                else:
                    division = self.division_3d
                meshbuilder = MeshBuilder(division)
                center = self.get_list3d(self.meshdict['Box']['Center'])
                mesh = meshbuilder.MeshBox(center, self.meshdict['Box']['Length'],
                                           self.meshdict['Box']['Width'],
                                           self.meshdict['Box']['Height'])

            if self.meshdict['Method'] == 'UnitCircleMesh':
                if 'Division' in self.meshdict:
                    division = self.get_list1d(self.meshdict['Division'])
                else:
                    division = self.division_1d
                meshbuilder = MeshBuilder(division)
                mesh = meshbuilder.MeshUnitCircle()

            if self.meshdict['Method'] == 'UnitSphereMesh':
                if 'Division' in self.meshdict:
                    division = self.get_list1d(self.meshdict['Division'])
                else:
                    division = self.division_1d
                meshbuilder = MeshBuilder(division)
                mesh = meshbuilder.MeshUnitSphere()

            if self.meshdict['Method'] == 'CircleMesh':
                if 'Division' in self.meshdict:
                    division = self.get_list1d(self.meshdict['Division'])
                else:
                    division = self.division_1d
                meshbuilder = MeshBuilder(division)
                mesh = meshbuilder.MeshCircle(self.meshdict['Circle']['Radius'])

            if self.meshdict['Method'] == 'SphereMesh':
                if 'Division' in self.meshdict:
                    division = self.get_list1d(self.meshdict['Division'])
                else:
                    division = self.division_1d
                meshbuilder = MeshBuilder(division)
                mesh = meshbuilder.MeshSphere(self.meshdict['Sphere']['Radius'])

        else:
            domain = self.domain_operation()
            if 'Division' in self.meshdict:
                division = self.get_list1d(self.meshdict['Division'])
            else:
                division = self.division_1d
            meshbuilder = MeshBuilder(division)
            mesh = meshbuilder.Generate(domain)

        return mesh

    def mesh_with_marker(self):
        key = list(self.meshdict['Domain'].items())[0][0]
        domain = self.domain_generator(key, self.meshdict['Domain'][key])
        if 'Subdomains' in self.meshdict:
            domain.set_subdomain(self.meshdict['Domain']['Marker'], domain)
            for sd in self.meshdict['Subdomains']:
                key = list(sd.items())[0][0]
                domain_sub = self.domain_generator(key, sd[key])
                domain.set_subdomain(sd['Marker'], domain_sub)

            if 'Division' in self.meshdict:
                division = self.get_list1d(self.meshdict['Division'])
            else:
                division = self.division_1d
            meshbuilder = MeshBuilder(division)
            mesh = meshbuilder.Generate(domain)
            d = mesh.topology().dim()  # dimensionality of the problem
            markers = MeshFunction("size_t", mesh, d, mesh.domains())

        return mesh, markers

    def create_mesh_object(self):
        if 'Input' in self.meshdict:
            meshinputer = MeshInputer(self.meshdict['Input'])
            mesh = meshinputer.mesh_input()
        else:
            if 'Domain' in self.meshdict:
                mesh, markers = self.mesh_with_marker()
                return mesh, markers
            else:
                mesh = self.mesh_builder()
        return mesh
