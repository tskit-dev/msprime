#
# Copyright (C) 2021 University of Oxford
#
# This file is part of msprime.
#
# msprime is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# msprime is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with msprime.  If not, see <http://www.gnu.org/licenses/>.
#
"""
Define formats used for simulation input as JSON and related formats.
"""
from __future__ import annotations

import copy
import dataclasses
import json

import demes
from ruamel.yaml import YAML

import msprime


@dataclasses.dataclass
class SimulationConfig:
    ancestry_kwargs: dict
    mutations_kwargs: dict = None


def parse_ancestry_json(data):
    data = copy.deepcopy(data)
    if "start_time" in data or "end_time" in data:
        raise ValueError(
            "specifying time values not currently supported as too confusing"
        )
    if "demography" in data:
        demes_dict = data["demography"]
        # TODO nasty going back to JSON here - can we make a demes.fromdict()
        # function to do this directly?
        demes_model = demes.loads(json.dumps(demes_dict), format="json")
        data["demography"] = msprime.Demography.from_demes(demes_model)
    return data


def parse_mutations_json(data):

    if "start_time" in data or "end_time" in data:
        raise ValueError(
            "specifying time values not currently supported as too confusing"
        )
    return data


def parse_yaml(text):

    yaml = YAML(typ="safe")
    data = yaml.load(text)
    config = SimulationConfig(parse_ancestry_json(data["ancestry"]))
    if "mutations" in data:
        config.mutations_kwargs = parse_mutations_json(data["mutations"])
    return config
