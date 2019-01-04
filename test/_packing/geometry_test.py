#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 Malcolm Ramsay <malramsay64@gmail.com>
#
# Distributed under terms of the MIT license.

import pytest

from _packing import Vect2, on_segment, segments_cross, triplet_orientation


@pytest.fixture
def points():
    """
               |     C
              1+     +
               |
               |
         A     |E    B
    -----+-----+-----+-----
        -1     |     1
               |
         D     |
         +   -1+
               |
    """
    from typing import NamedTuple

    class Vects(NamedTuple):
        a: Vect2 = Vect2(-1, 0)
        b: Vect2 = Vect2(1, 0)
        c: Vect2 = Vect2(1, 1)
        d: Vect2 = Vect2(-1, -1)
        e: Vect2 = Vect2(0, 0)

    return Vects()


def test_triplet_orientation(points):
    # Clockwise
    assert triplet_orientation(points.c, points.b, points.a) == 1
    assert triplet_orientation(points.c, points.d, points.a) == 1
    assert triplet_orientation(points.b, points.d, points.a) == 1

    # Anticlockwise
    assert triplet_orientation(points.a, points.b, points.c) == -1
    assert triplet_orientation(points.e, points.b, points.c) == -1
    assert triplet_orientation(points.d, points.b, points.e) == -1

    # Colinear
    assert triplet_orientation(points.a, points.e, points.b) == 0
    assert triplet_orientation(points.d, points.e, points.c) == 0


def test_on_segment(points):
    assert on_segment(points.a, points.e, points.b) is True
    assert on_segment(points.a, points.b, points.e) is False
    assert on_segment(points.d, points.e, points.c) is True
    assert on_segment(points.e, points.d, points.c) is False


def test_segments_cross(points):
    assert segments_cross(points.a, points.b, points.c, points.d) is True
    assert segments_cross(points.e, points.b, points.c, points.d) is True
    assert segments_cross(points.a, points.d, points.b, points.c) is False
