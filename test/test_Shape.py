#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 Malcolm Ramsay <malramsay64@gmail.com>
#
# Distributed under terms of the MIT license.

""""""

import pytest

from _packing import Shape


def test_init():
    Shape("test", [1, 1, 1, 1], 0, 0)


@pytest.fixture
def shape():
    return Shape("test", [1, 1, 1, 1], 0, 0)


def test_get_points(shape):
    assert shape.get_point(0) == 1
    assert shape.get_point(1) == 1
    assert shape.get_point(2) == 1
    assert shape.get_point(3) == 1
    with pytest.raises(IndexError):
        shape.get_point(4)
    with pytest.raises(IndexError):
        shape.get_point(-1)
