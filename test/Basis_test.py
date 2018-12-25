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

from _packing import Basis, CellAngleBasis, CellLengthBasis, FixedBasis


class BasisFixture(NamedTuple):
    expected_value: float
    expected_min_value: float
    expected_max_value: float
    basis: Basis
    cell_x: Optional[Basis] = None
    cell_y: Optional[Basis] = None


@pytest.fixture(params=["Basis", "CellAngleBasis", "CellLengthBasis"])
def basis_fixture(request):
    value = 0.5
    min_val = 0
    max_val = 1
    cell_x = None
    cell_y = None
    if request.param == "Basis":
        basis = Basis(value, min_val, max_val)
    elif request.param == "CellAngleBasis":
        step = 0.01
        cell_x = Basis(0.5, 0, 1)
        cell_y = Basis(0.5, 0, 1)
        basis = CellAngleBasis(value, min_val, max_val, step, cell_x, cell_y)
    elif request.param == "FixedBasis":
        min_val = value
        max_val = value
        basis = FixedBasis(value)
    elif request.param == "CellLengthBasis":
        basis = CellLengthBasis(value, min_val, max_val)

    else:
        raise NotImplementedError(f"The basis {request.param} is not implemented.")

    return BasisFixture(
        expected_value=value,
        expected_min_value=min_val,
        expected_max_value=max_val,
        basis=basis,
        cell_x=cell_x,
        cell_y=cell_y,
    )


def test_init(basis_fixture):
    basis = basis_fixture.basis
    assert isinstance(basis, Basis)


def test_value(basis_fixture):
    basis = basis_fixture.basis
    assert basis.value == basis_fixture.expected_value


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


@pytest.mark.parametrize("basis_fixture", ["FixedBasis"], indirect=["basis_fixture"])
def test_fixed_basis(basis_fixture):
    basis = basis_fixture.basis
    assert basis.value == basis_fixture.expected_value

    basis.value = 0
    assert basis.value == basis_fixture.expected_value
    assert basis.value_range() == 0


@pytest.mark.parametrize(
    "basis_fixture", ["CellAngleBasis"], indirect=["basis_fixture"]
)
def test_cell_angle_basis(basis_fixture):
    basis = basis_fixture.basis
    assert basis_fixture.cell_x is basis.cell_x_len
    assert basis_fixture.cell_y is basis.cell_y_len

    # Make angle smaller
    basis.value = 0.4
    assert basis_fixture.cell_x.value > 0.5
    assert basis_fixture.cell_y.value > 0.5

    # Reset angle
    basis.reset_value()
    assert math.isclose(basis.value, basis_fixture.expected_value)
    assert basis_fixture.cell_x.value == 0.5
    assert basis_fixture.cell_y.value == 0.5

    # Make angle larger
    basis.value = 0.6
    assert basis_fixture.cell_x.value < 0.5
    assert basis_fixture.cell_y.value < 0.5
