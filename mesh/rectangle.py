import fenics
import dolfin
import mshr
class Rectangle:
    def __init__(self, inputs):
        pt1 = inputs['End Points'][0]
        pt2 = inputs['End Points'][1]
        self.object = mshr.Rectangle(fenics.Point(pt1),fenics.Point(pt2))

