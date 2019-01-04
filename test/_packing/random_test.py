#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 Malcolm Ramsay <malramsay64@gmail.com>
#
# Distributed under terms of the MIT license.

from _packing import fluke


def test_fluke():
    assert isinstance(fluke(), float)

def test_range():
    for i in range(10000):
        val = fluke()
        assert val < 1
        assert val > 0

