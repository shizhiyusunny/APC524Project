#!/usr/bin/env python
# coding: utf-8

from ..material import material
import pytest


def test_material1():
    material_constants = {'E': 0.5, 'nu': 0.3}
    material.MaterialConstant.processor(material_constants)
    assert abs(material_constants['mu']-0.1923) < 0.0001


def test_material2():
    material_constants = {'E': 0.6, 'nu': 0.2}
    material.MaterialConstant.processor(material_constants) 
    assert abs(material_constants['lambda']-0.125) < 0.0001

