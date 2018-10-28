#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 Malcolm Ramsay <malramsay64@gmail.com>
#
# Distributed under terms of the MIT license.

import os
import sys

from distutils.core import Extension
from distutils.command.build_ext

import numpy
from cython import Cythonize

extensions = []
