#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 Malcolm Ramsay <malramsay64@gmail.com>
#
# Distributed under terms of the MIT license.

"""

"""

import re

import numexpr
from ruamel.yaml import YAML, yaml_object

yaml = YAML(typ="unsafe")


@yaml_object(yaml)
class WyckoffInfo:
    yaml_tag: str = "WyckoffInfo"
    multiplicity: int
    letter: str
    rot: int
    mirror_x: bool
    mirror_y: bool
    coordinates: str

    def __init__(self, multiplicity, letter, symmetries, coordinates):
        self.multiplicity = int(multiplicity)
        self.letter = str(letter)
        self.rot, self.mirror_x, self.mirror_y = self.parse_symmetries(symmetries)
        self.coordinates = self.parse_coordinates(coordinates)

    @staticmethod
    def parse_symmetries(symmetry):
        if "m" not in (symmetry) and "." not in symmetry:
            return (int(symmetry), False, False)
        symmetry = list(symmetry)
        if symmetry[0] == ".":
            symmetry[0] = 1
        return (int(symmetry[0]), "m" in symmetry[1], "m" in symmetry[2])

    @staticmethod
    def parse_coordinates(coordinates):
        coord_list = re.findall('\[[^\]]*\]|\([^\)]*\)|"[^"]*"|\S+', coordinates)
        return [SymmetryOps(cmd) for cmd in coord_list]

    def __repr__(self):
        return (
            f"WyckoffInfo(letter='{self.letter}', "
            f"mult={self.multiplicity},  rot={self.rot}, "
            f"mirror_x={self.mirror_x}, mirror_y={self.mirror_y}, "
            f"coordinates='{self.coordinates}')"
        )

    @classmethod
    def from_row(cls, row):
        return WyckoffInfo(*[element.string for element in row.find_all("td")])


@yaml_object(yaml)
class SymmetryOps:
    def __init__(self, string):
        self.x_op, self.y_op = string.strip("()").split(",")
        self.x_op = self.x_op.strip()
        self.y_op = self.y_op.strip()

    def convert(self, x, y):
        return numexpr.evaluate(self.x_op), numexpr.evaluate(self.y_op)

    def __repr__(self):
        return f"SymOp(x={self.x_op}, y={self.y_op})"

    def create_matrix(self):
        c = self.convert(0, 0)
        x = self.convert(1, 0)
        y = self.convert(0, 1)
        x_vals = [x[0] - c[0], y[0] - c[0], c[0]]
        y_vals = [x[1] - c[1], y[1] - c[1], c[1]]
        return [x_vals, y_vals]

    def matrix_convert(self, x, y):
        mat = self.create_matrix()
        x_ret = mat[0][0] * x + mat[0][1] * y + mat[0][2]
        y_ret = mat[1][0] * x + mat[1][1] * y + mat[1][2]
        return (x_ret, y_ret)
