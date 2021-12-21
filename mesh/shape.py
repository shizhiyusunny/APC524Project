from fenics import *
from mshr import *


class Shape:

    def __init__(self, center):
        self.d = len(center)
        if self.d == 1:
            self.centerx = center[0]
        if self.d == 2:
            self.centerx = center[0]
            self.centery = center[1]
        if self.d == 3:
            self.centerx = center[0]
            self.centery = center[1]
            self.centerz = center[2]

    # 1D shape
    def interval(self, length):
        pt1x = self.centerx - 0.5 * length
        pt2x = self.centerx + 0.5 * length
        domain = Interval(pt1x, pt2x)
        return domain

    # 2D shape
    def circle(self, r):
        domain = Circle(Point(self.centerx, self.centery), r)
        return domain

    def rectangle(self, length, breadth):
        pt1x = self.centerx - 0.5 * length
        pt1y = self.centery - 0.5 * breadth
        pt2x = self.centerx + 0.5 * length
        pt2y = self.centery + 0.5 * breadth
        domain = Rectangle(Point(pt1x, pt1y), Point(pt2x, pt2y))
        return domain

    def square(self, length):
        domain = self.rectangle(length, length)
        return domain

    # 3D shape
    def sphere(self, r):
        domain = Sphere(Point(self.centerx, self.centery, self.centerz), r)
        return domain

    def box(self, length, width, height):
        pt1x = self.centerx - 0.5 * length
        pt1y = self.centery - 0.5 * width
        pt1z = self.centerz - 0.5 * height
        pt2x = self.centerx + 0.5 * length
        pt2y = self.centery + 0.5 * width
        pt2z = self.centerz + 0.5 * height
        domain = Box(Point(pt1x, pt1y, pt1z), Point(pt2x, pt2y, pt2z))
        return domain

    def cube(self, length):
        domain = self.box(length, length, length)
        return domain
