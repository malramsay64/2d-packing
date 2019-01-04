#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 Malcolm Ramsay <malramsay64@gmail.com>
#
# Distributed under terms of the MIT license.

"""Define a series of shapes for the packing algorithm.

All shapes are derived from the py:class:`Shape` base class, 
and additional shapes can be created by sub-classing.

"""
from abc import ABC
from typing import List, Optional

from _packing import Shape


class Polygon(Shape):
    sides: int

    def __init__(self, sides, name=None):
        if sides < 3:
            raise ValueError(f"A polygon must have 3 or more sides, got {sides}.")
        self.sides = sides
        if name is None:
            name = f"Polygon-{sides}"

        rotational_symmetries = self.sides
        mirrors = sides if sides % 2 == 0 else self.sides // 2
        super().__init__(name, [1] * self.sides, rotational_symmetries, mirrors)


class Square(Polygon):
    def __init__(self):
        super().__init__(4, "Square")


class Triangle(Polygon):
    def __init__(self):
        super().__init__(3, "Triangle")


class Circle(Polygon):
    def __init__(self, resolution=360):
        super().__init__(resolution, "Circle")
