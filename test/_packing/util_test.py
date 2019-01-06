#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 Malcolm Ramsay <malramsay64@gmail.com>
#
# Distributed under terms of the MIT license.

import itertools
from collections import OrderedDict

import pytest
from hypothesis import given
from hypothesis.strategies import floats, integers, lists

from _packing import combinations, uniqueify

INT_MAX = 2 ^ 63 - 1


def unique_combs(values):
    vals = OrderedDict.fromkeys([tuple(i) for i in values])
    return [list(i) for i in vals.keys()]


def reference_unique(values):
    vals = OrderedDict.fromkeys([tuple(i) for i in values])
    return [list(i) for i in vals.keys()]


def reference_combinations(values, take):
    return [list(val) for val in itertools.combinations(values, take)]


@given(lists(integers(max_value=INT_MAX, min_value=-INT_MAX), min_size=3, unique=True))
def test_combinations_int_unique(values):
    take = 3
    reference_comb = reference_combinations(values, take)
    comb = combinations(values, take)
    assert comb == reference_comb


@given(lists(integers(max_value=INT_MAX, min_value=-INT_MAX), min_size=3, unique=False))
def test_combinations_int(values):
    take = 3
    reference_comb = reference_combinations(values, take)
    comb = combinations(values, take)
    assert comb == reference_comb


@given(lists(floats(allow_nan=False), min_size=3, unique=True))
def test_combinations_double_unique(values):
    take = 3
    reference_comb = reference_combinations(values, take)
    comb = combinations(values, take)
    assert comb == reference_comb


@given(lists(floats(allow_nan=False), min_size=3, unique=False))
def test_combinations_double(values):
    take = 3
    reference_comb = reference_combinations(values, take)
    comb = combinations(values, take)
    assert comb == reference_comb


@given(lists(lists(floats(allow_nan=False), min_size=1, unique=False), min_size=1))
def test_uniqueify_float(values):
    assert uniqueify(values) == reference_unique(values)


@given(
    lists(
        lists(
            integers(max_value=INT_MAX, min_value=-INT_MAX), min_size=1, unique=False
        ),
        min_size=1,
    )
)
def test_uniqueify_int(values):
    assert uniqueify(values) == reference_unique(values)
