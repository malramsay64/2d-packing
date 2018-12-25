#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 Malcolm Ramsay <malramsay64@gmail.com>
#
# Distributed under terms of the MIT license.

from typing import NamedTuple

import pytest
from hypothesis import given
from hypothesis.strategies import floats

from _packing import Basis


class BasisFixture(NamedTuple):
    expected_value: float
    expected_min_value: float
    expected_max_value: float
    basis: Basis


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
