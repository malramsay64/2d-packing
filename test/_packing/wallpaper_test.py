#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 Malcolm Ramsay <malramsay64@gmail.com>
#
# Distributed under terms of the MIT license.

import pytest

from _packing import Mirror, SymmetryTransform, Vect2, Vect3


def mirror_states():
    return [func for func in dir(Mirror) if not func.startswith("_")]


@pytest.mark.parametrize("state", mirror_states())
def test_Mirror(state):
    angle = int(state[1:])
    assert Mirror(angle) == getattr(Mirror, state)


def test_SymmetryTransform():
    st = SymmetryTransform(Vect3(1, 0, 0), Vect3(0, 1, 0), 0, Mirror(0))
    assert st.real_to_fractional(Vect3(0.5, 0.5, 0)) == Vect2(0.5, 0.5)
