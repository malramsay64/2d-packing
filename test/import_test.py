#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 Malcolm Ramsay <malramsay64@gmail.com>
#
# Distributed under terms of the MIT license.

"""Ensure the importing of packages is working correctly."""


def test_packing():
    import _packing

def test_pypacking():
    import pypacking

def test_pypacking_packing():
    from pypacking import _packing
