#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 Malcolm Ramsay <malramsay64@gmail.com>
#
# Distributed under terms of the MIT license.

""""""

import pytest

from pypacking import shapes

all_shapes = {
    "Shape": shapes.Shape("test", [1] * 10, 0, 0),
    "Circle": shapes.Circle(),
    "Triangle": shapes.Triangle(),
    "Square": shapes.Square(),
}

for i in range(3, 13):
    all_shapes[f"Polygon-{i}"] = shapes.Polygon(i)


@pytest.fixture(params=all_shapes.values(), ids=all_shapes.keys())
def shape(request):
    yield request.param


def test_init(shape):
    assert isinstance(shape, shapes.Shape)


@pytest.mark.parametrize("sides", range(3))
def test_polygon_values_fail(sides):
    with pytest.raises(ValueError):
        shapes.Polygon(sides)
