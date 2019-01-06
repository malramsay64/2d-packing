#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 Malcolm Ramsay <malramsay64@gmail.com>
#
# Distributed under terms of the MIT license.

import itertools

import pytest
from hypothesis import given
from hypothesis.strategies import floats, integers, lists

from _packing import py_combinations

INT_MAX = 2 ^ 63 - 1


@given(lists(integers(max_value=INT_MAX, min_value=-INT_MAX), min_size=3, unique=True))
def test_combinations_int(values):
    take = 3
    iter_comb = [list(val) for val in itertools.combinations(values, take)]
    py_comb = py_combinations(values, take)
    assert py_comb == iter_comb


@given(lists(floats(allow_nan=False), min_size=3, unique=True))
def test_combinations_double(values):
    take = 3
    iter_comb = [list(val) for val in itertools.combinations(values, take)]
    py_comb = py_combinations(values, take)
    assert py_comb == iter_comb
