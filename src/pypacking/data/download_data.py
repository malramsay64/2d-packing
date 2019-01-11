#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 Malcolm Ramsay <malramsay64@gmail.com>
#
# Distributed under terms of the MIT license.

"""This file downloads the data for all the 2D space groups.

The bilbao crystallographic server contains information on all the different wallpaper
groups and space groups, however that data is not available in an easy to read or
complete format. This module downloads the data for all the Wyckoff sites for each
wallpaper group, getting the data into a format that is easy to convert between
representations. Part of this conversion is the creation of instances of
py:class:`WyckoffInfo` and py:class:`SymmetryOps` which can be dumped and loaded from a
file, and have a number of transformations implemented.

"""

import re
from pathlib import Path

import requests
from bs4 import BeautifulSoup

from .wyckoff import SymmetryOps, WyckoffInfo, yaml


def parse_page(content):
    soup = BeautifulSoup(content, "html.parser")
    re_string = r"^Wyckoff Positions of Plane Group (?P<group>.*) \(.*\)$"
    wallpaper_group = re.sub(re_string, "\g<group>", soup.h2.text).replace(" ", "")
    table = soup.find(lambda tag: tag.name == "table" and tag.has_attr("border"))
    return [WyckoffInfo.from_row(wallpaper_group, row) for row in table.find_all("tr")]


def get_wallpaper_groups():
    url = "http://www.cryst.ehu.es/cgi-bin/plane/programs/nph-plane_wp-list"
    sites = []
    for i in range(1, 18):
        req = requests.put(url, data={"gnum": str(i)})
        sites.extend(parse_page(req.content))
    return sites


def save_wallpaper_groups():
    sites = get_wallpaper_groups()

    with open(Path(__file__).parent / "wallpaper_groups.yaml", "w") as dst:
        yaml.dump(sites, dst)


def load_wallpaper_groups():
    with open(Path(__file__).parent / "wallpaper_groups.yaml", "r") as src:
        return yaml.load(src)


if __name__ == "__main__":
    save_wallpaper_groups()
