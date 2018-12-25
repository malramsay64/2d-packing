#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 Malcolm Ramsay <malramsay64@gmail.com>
#
# Distributed under terms of the MIT license.

""""""

import math

import pytest

from _packing import Shape


def test_init():
    Shape("test", [1, 1, 1, 1], 0, 0)


def get_shape_name(sides: int):
    return {
        3: "Triangle",
        4: "Square",
        5: "Pentagon",
        6: "Hexagon",
        7: "Heptagon",
        8: "Octagon",
        9: "Nonagon",
        10: "Decagon",
    }.get(sides, "Polygon")


@pytest.fixture(params=range(3, 11), ids=[get_shape_name(i) for i in range(3, 11)])
def polygon(request):
    sides = request.param
    return sides, Shape(get_shape_name(sides), [1] * sides, 0, 0)


def test_get_points(polygon):
    sides, shape = polygon
    for i in range(sides):
        assert shape.get_point(i) == 1

    with pytest.raises(IndexError):
        shape.get_point(-1)
    with pytest.raises(IndexError):
        shape.get_point(sides)


def test_resolution(polygon):
    sides, shape = polygon
    assert shape.resolution() == sides


def test_angular_step(polygon):
    sides, shape = polygon
    assert math.isclose(shape.angular_step(), math.tau / sides)


def test_area(polygon):
    sides, shape = polygon
    area = 0.5 * math.sin(math.tau / sides) * sides
    assert math.isclose(shape.area(), area)
