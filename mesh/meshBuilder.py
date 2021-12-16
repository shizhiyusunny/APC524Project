from fenics import *
import mshr
from . import shape
class MeshBuilder:
    def build(inputs):
        if 'Domain' in inputs:
            domain = shape.Shape.build(inputs['Domain'])
            if 'Subdomains' in inputs:
                domain.object.set_subdomain(inputs['Domain']['Marker'], domain.object)
                for sd in inputs['Subdomains']:
                    domain.object.set_subdomain(sd['Marker'], shape.Shape.build(sd).object)
            if 'divisions' not in inputs:
                inputs['divisions'] = 100
            mesh = mshr.generate_mesh(domain.object, inputs['divisions'])
            d = mesh.topology().dim() # dimensionality of the problem
            markers = MeshFunction("size_t", mesh, d , mesh.domains())
        return mesh, markers
