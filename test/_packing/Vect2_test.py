#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 Malcolm Ramsay <malramsay64@gmail.com>
#
# Distributed under terms of the MIT license.

import math
import operator

import pytest
from hypothesis import given
from hypothesis.strategies import floats

from _packing import Vect2


def test_init():
    assert Vect2()
    assert Vect2(1, 1)


def test_access():
    v = Vect2(0, 1)
    assert v.x == 0
    assert v.y == 1


@pytest.mark.parametrize("operator", [operator.add, operator.mul, operator.sub])
@given(floats(), floats())
def test_operators(operator, x, y):
    v1 = Vect2(x, y)
    v2 = Vect2(y, x)
    v_op = operator(v1, v2)
    if math.isnan(operator(x, y)):
        assert math.isnan(v_op.x)
        assert math.isnan(v_op.y)
    else:
        assert v_op.x == operator(x, y)
        assert v_op.y == operator(y, x)


@given(floats(), floats())
def test_norm(x, y):
    v = Vect2(x, y)
    if math.isnan(x) or math.isnan(y):
        assert math.isnan(v.norm())
    else:
        assert v.norm() == math.sqrt(x * x + y * y)


@given(floats(), floats())
def test_norm_sq(x, y):
    v = Vect2(x, y)
    if math.isnan(x) or math.isnan(y):
        assert math.isnan(v.norm_sq())
    else:
        assert v.norm_sq() == x * x + y * y
