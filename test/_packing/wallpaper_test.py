#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 Malcolm Ramsay <malramsay64@gmail.com>
#
# Distributed under terms of the MIT license.

import pytest

from _packing import Mirror


def mirror_states():
    return [func for func in dir(Mirror) if not func.startswith("_")]


@pytest.mark.parametrize("state", mirror_states())
def test_Mirror(state):
    angle = int(state[1:])
    assert Mirror(angle) == getattr(Mirror, state)
