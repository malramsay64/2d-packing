#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 Malcolm Ramsay <malramsay64@gmail.com>
#
# Distributed under terms of the MIT license.

import pytest

from _packing import (
    Mirror,
    SymmetryTransform,
    Vect2,
    Vect3,
    WallpaperGroup,
    WyckoffSite,
)


def mirror_states():
    return [func for func in dir(Mirror) if not func.startswith("_")]


@pytest.mark.parametrize("state", mirror_states())
def test_Mirror(state):
    angle = int(state[1:])
    assert Mirror(angle) == getattr(Mirror, state)


def test_SymmetryTransform():
    st = SymmetryTransform(Vect3(1, 0, 0), Vect3(0, 1, 0))
    assert st.real_to_fractional(Vect2(0.5, 0.5)) == Vect2(0.5, 0.5)


@pytest.fixture(params=[-1, 1])
def symmetry(request):
    sign = request.param
    yield SymmetryTransform(Vect3(sign * 1, 0, 0), Vect3(0, sign * 1, 0))


def test_WyckoffSite(symmetry):
    site = WyckoffSite("a", [symmetry])
    assert site.multiplicity() == 1
    assert site.vary_x() is True
    assert site.vary_x() is True
    assert site.mirror_type() == 0


@pytest.fixture
def wyckoff_site(symmetry):
    yield WyckoffSite("a", [symmetry])


def test_WallpaperGroup(wyckoff_site):
    wallpaper = WallpaperGroup("p1", [wyckoff_site])
    assert wallpaper.label == "p1"
