#!/usr/bin/env python
# coding: utf-8


from ..mesh import shape, setshape
from fenics import *
from mshr import *
import pytest


def test_shape3d():
    center = [0,0,0]
    shape1 = shape.Shape(center)
    domain1 = shape1.box(3,2,4)
    shape2 = shape.Shape(center)
    domain2 = shape2.sphere(0.1)
    sets = setshape.SetShape(domain1,domain2)
    domain = sets.subtraction() 
    domain3 = Box(Point(-1.5,-1,-2),Point(1.5,1,2)) - Sphere(Point(0,0,0),0.5)
    assert type(domain) == type(domain3)


def test_shape2d():
    center = [0,0]
    shape1 = shape.Shape(center)
    domain1 = shape1.rectangle(2,3)
    shape2 = shape.Shape(center)
    domain2 = shape2.circle(0.5)
    sets = setshape.SetShape(domain1,domain2)
    domain = sets.subtraction()
    domain3 = Rectangle(Point(-1,-1.5),Point(1,1.5)) - Circle(Point(0,0),0.5)
    assert type(domain) == type(domain3)
