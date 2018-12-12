#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 Malcolm Ramsay <malramsay64@gmail.com>
#
# Distributed under terms of the MIT license.

"""

"""
from math import pi
from typing import List

import attr


@attr.s(auto_attribs=True)
class Image:
    x_coeffs: List[float]
    y_coeffs: List[float]

    rotation_offset: float = 0
    flipped: bool = False

    site_mirror_0: bool = False
    site_mirror_90: bool = False
    site_mirror_45: bool = False
    site_mirror_135: bool = False
    site_mirror_30: bool = False
    site_mirror_60: bool = False
    site_mirror_330: bool = False
    site_mirror_300: bool = False


@attr.s(auto_attribs=True)
class WyckoffType:
    letter: str
    site_rotations: int
    site_mirrors: int
    images: List[Image]
    variability: int = 0

    @property
    def multiplicity(self):
        return len(self.images)


@attr.s(auto_attribs=True)
class WallpaperGroup:
    label: str
    wyckoffs: List[WyckoffType]
    num_symmetries: int
    a_b_equal: bool = False
    hexagonal: bool = False
    rectangular: bool = False

    @property
    def num_wyckoffs(self):
        return len(self.wyckoffs)


wallpaper_p1 = WallpaperGroup(
    label="p1",
    num_symmetries=1,
    wyckoffs=[
        WyckoffType(
            letter="a",
            site_rotations=1,
            site_mirrors=1,
            images=[Image(x_coeffs=[1.0, 0.0, 0.0], y_coeffs=[0.0, 1.0, 0.0])],
        )
    ],
)

wallpaper_p2 = WallpaperGroup(
    label="p2",
    num_symmetries=2,
    wyckoffs=[
        WyckoffType(
            letter="a",
            site_rotations=2,
            site_mirrors=0,
            images=[Image(x_coeffs=[0.0, 0.0, 0.0], y_coeffs=[0.0, 0.0, 0.0])],
        ),
        WyckoffType(
            letter="b",
            site_rotations=2,
            site_mirrors=0,
            images=[Image(x_coeffs=[0.0, 0.0, 0.0], y_coeffs=[0.0, 0.0, 0.5])],
        ),
        WyckoffType(
            letter="c",
            site_rotations=2,
            site_mirrors=0,
            images=[Image(x_coeffs=[0.0, 0.0, 0.5], y_coeffs=[0.0, 0.0, 0.0])],
        ),
        WyckoffType(
            letter="d",
            site_rotations=2,
            site_mirrors=0,
            images=[Image(x_coeffs=[0.0, 0.0, 0.5], y_coeffs=[0.0, 0.0, 0.5])],
        ),
        WyckoffType(
            letter="e",
            site_rotations=1,
            site_mirrors=0,
            variability=1,
            images=[
                Image(
                    x_coeffs=[1.0, 0.0, 0.0],
                    y_coeffs=[0.0, 1.0, 0.0],
                    rotation_offset=0,
                ),
                Image(
                    x_coeffs=[-1.0, 0.0, 0.0],
                    y_coeffs=[0.0, -1.0, 0.0],
                    rotation_offset=pi,
                ),
            ],
        ),
    ],
)

wallpaper_pm = WallpaperGroup(
    label="pm",
    rectangular=True,
    num_symmetries=2,
    wyckoffs=[
        WyckoffType(
            letter="a",
            variability=1,
            site_rotations=1,
            site_mirrors=1,
            images=[
                Image(
                    x_coeffs=[0.0, 0.0, 0.0],
                    y_coeffs=[0.0, 1.0, 0.0],
                    site_mirror_90=True,
                )
            ],
        ),
        WyckoffType(
            letter="b",
            variability=1,
            site_rotations=1,
            site_mirrors=1,
            images=[
                Image(
                    x_coeffs=[0.0, 0.0, 0.5],
                    y_coeffs=[0.0, 1.0, 0.0],
                    site_mirror_90=True,
                )
            ],
        ),
        WyckoffType(
            letter="c",
            variability=1,
            site_rotations=1,
            images=[
                Image(
                    x_coeffs=[1.0, 0.0, 0.0],
                    y_coeffs=[0.0, 1.0, 0.0],
                    rotation_offset=0,
                    flipped=False,
                ),
                Image(
                    x_coeffs=[-1.0, 0.0, 0.0],
                    y_coeffs=[0.0, 1.0, 0.0],
                    rotation_offset=pi,
                    flipped=True,
                ),
            ],
        ),
    ],
)

wallpaper_pg = WallpaperGroup(
    label="pg",
    rectangular=True,
    num_symmetries=2,
    wyckoffs=[
        WyckoffType(
            letter="a",
            variability=1,
            site_rotations=1,
            images=[
                Image(
                    x_coeffs=[1.0, 0.0, 0.0],
                    y_coeffs=[0.0, 1.0, 0.0],
                    rotation_offset=0,
                    flipped=False,
                ),
                Image(
                    x_coeffs=[-1.0, 0.0, 0.0],
                    y_coeffs=[0.0, 1.0, 0.5],
                    rotation_offset=pi,
                    flipped=True,
                ),
            ],
        )
    ],
)

wallpaper_cm = WallpaperGroup(
    label="cm",
    rectangular=True,
    num_symmetries=4,
    wyckoffs=[
        WyckoffType(
            letter="a",
            variability=True,
            site_rotations=1,
            site_mirrors=1,
            images=[
                Image(
                    x_coeffs=[0.0, 0.0, 0.0],
                    y_coeffs=[0.0, 1.0, 0.0],
                    site_mirror_90=True,
                ),
                Image(
                    x_coeffs=[0.0, 0.0, 0.5],
                    y_coeffs=[0.0, 1.0, 0.5],
                    site_mirror_90=True,
                ),
            ],
        ),
        WyckoffType(
            letter="b",
            variability=True,
            site_rotations=1,
            site_mirrors=0,
            images=[
                Image(
                    x_coeffs=[1.0, 0.0, 0.0],
                    y_coeffs=[0.0, 1.0, 0.0],
                    rotation_offset=0,
                    flipped=False,
                ),
                Image(
                    x_coeffs=[-1.0, 0.0, 0.0],
                    y_coeffs=[0.0, 1.0, 0.0],
                    rotation_offset=pi,
                    flipped=True,
                ),
                Image(
                    x_coeffs=[1.0, 0.0, 0.5],
                    y_coeffs=[0.0, 1.0, 0.5],
                    rotation_offset=0,
                    flipped=False,
                ),
                Image(
                    x_coeffs=[-1.0, 0.0, 0.5],
                    y_coeffs=[0.0, 1.0, 0.5],
                    rotation_offset=pi,
                    flipped=True,
                ),
            ],
        ),
    ],
)

wallpaper_groups = [
    wallpaper_p1,
    wallpaper_p2,
    wallpaper_pg,
    wallpaper_pm,
    wallpaper_cm,
]
