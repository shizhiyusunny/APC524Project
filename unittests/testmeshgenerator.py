#!/usr/bin/env python
# coding: utf-8


from fenics import *
from mshr import *
from ..mesh import meshmenerator
import pytest


def test_meshgenerator():
    meshdict = {'Rectangle':{'Center':'0.0, 0.0', 'Length': 1.0, 'Breadth': 1.0}}
    meshg = meshmenerator.MeshGenerator(meshdict)
    mesh1 = meshg.create_mesh_object()
    domain = Rectangle(Point(-0.5,-0.5),Point(0.5,0.5))
    mesh2 = generate_mesh(domain,20)
    assert type(mesh1) == type(mesh2)

