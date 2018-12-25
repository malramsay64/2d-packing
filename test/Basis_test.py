#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 Malcolm Ramsay <malramsay64@gmail.com>
#
# Distributed under terms of the MIT license.

from typing import NamedTuple

import pytest

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
