#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 Malcolm Ramsay <malramsay64@gmail.com>
#
# Distributed under terms of the MIT license.

import math
from typing import NamedTuple, Optional

import pytest
from hypothesis import given
from hypothesis.strategies import floats

from _packing import Basis, CellAngleBasis, FixedBasis


class BasisFixture(NamedTuple):
    expected_value: float
    expected_min_value: float
    expected_max_value: float
    basis: Basis
    cell_x: Optional[Basis] = None
    cell_y: Optional[Basis] = None


@pytest.fixture(params=range(10))
def basis_fixture(request):
    yield BasisFixture(
        expected_value=request.param / 10,
        expected_min_value=0,
        expected_max_value=1,
        basis=Basis(request.param / 10, 0, 1),
    )


def test_init(basis_fixture):
    assert isinstance(basis_fixture.basis, Basis)


def test_value(basis_fixture):
    assert basis_fixture.basis.value == basis_fixture.expected_value


@given(floats(min_value=0, max_value=1))
def test_set_value(basis_fixture, new_value):
    basis = basis_fixture.basis
    basis.value = new_value
    assert basis.value == new_value


@given(floats(allow_nan=False))
def test_set_value_any(basis_fixture, new_value):
    basis = basis_fixture.basis
    basis.value = new_value
    if new_value < basis_fixture.expected_min_value:
        assert basis.value == basis_fixture.expected_min_value
    elif new_value > basis_fixture.expected_max_value:
        assert basis.value == basis_fixture.expected_max_value
    else:
        assert basis.value == new_value


@given(floats(allow_nan=False))
def test_reset_value_any(basis_fixture, new_value):
    basis = basis_fixture.basis
    basis.value = new_value
    basis.reset_value()
    assert basis.value == basis_fixture.expected_value


@pytest.fixture()
def fixed_basis_fixture():
    yield BasisFixture(
        expected_value=1,
        expected_min_value=1,
        expected_max_value=1,
        basis=FixedBasis(1),
    )


def test_fixed_basis(fixed_basis_fixture):
    basis = fixed_basis_fixture.basis
    assert basis.value == fixed_basis_fixture.expected_value

    basis.value = 0
    assert basis.value == fixed_basis_fixture.expected_value
    assert basis.value_range() == 0


@pytest.fixture()
def cell_angle_basis_fixture():
    value = math.pi / 3
    max_val = 3 * math.pi / 4
    min_val = math.pi / 4
    step = 0.01
    cell_x = Basis(0.5, 0, 1)
    cell_y = Basis(0.5, 0, 1)
    basis = BasisFixture(
        expected_value=value,
        expected_min_value=min_val,
        expected_max_value=max_val,
        basis=CellAngleBasis(value, min_val, max_val, step, cell_x, cell_y),
        cell_x=cell_x,
        cell_y=cell_y,
    )
    yield basis


def test_cell_angle_basis(cell_angle_basis_fixture):
    basis = cell_angle_basis_fixture.basis
    assert cell_angle_basis_fixture.cell_x is basis.cell_x_len
    assert cell_angle_basis_fixture.cell_y is basis.cell_y_len

    basis.value = math.pi / 4
    assert cell_angle_basis_fixture.cell_x.value > 0.5
    assert cell_angle_basis_fixture.cell_y.value > 0.5

    basis.reset_value()
    assert math.isclose(basis.value, cell_angle_basis_fixture.expected_value)
    assert cell_angle_basis_fixture.cell_x.value == 0.5
    assert cell_angle_basis_fixture.cell_y.value == 0.5

    basis.value = math.pi / 2
    assert cell_angle_basis_fixture.cell_x.value < 0.5
    assert cell_angle_basis_fixture.cell_y.value < 0.5
