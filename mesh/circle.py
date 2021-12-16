import mshr
import fenics

class Circle:
    def __init__(self, inputs):
        center = fenics.Point(inputs['Center'])
        radius = inputs['Radius']
        self.object = mshr.Circle(center, radius)
