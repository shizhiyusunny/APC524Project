#!/usr/bin/env python
# coding: utf-8


import pytest
from ..mesh import meshbuilder
from fenics import *
from mshr import *


def test_builder1d():
    division = [20]
    meshb = meshbuilder.MeshBuilder(division)
    mesh1 = meshb.UnitHyperCube()
    mesh2 = UnitIntervalMesh(*division)
    assert type(mesh1) == type(mesh2)


def test_builder2d():
    division = [20, 20]
    center = [0, 0]
    meshb = meshbuilder.MeshBuilder(division)
    mesh1 = meshb.MeshRect(center,3,2)
    mesh2 = RectangleMesh(Point(-1.5, -1), Point(1.5, 1), *division)
    assert type(mesh1) == type(mesh2)

