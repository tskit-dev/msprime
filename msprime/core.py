#
# Copyright (C) 2015-2021 University of Oxford
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
Core functions and classes used throughout msprime.
"""
from __future__ import annotations

import dataclasses
import numbers
import os
import random
import textwrap
from typing import Any
from typing import Dict
from typing import List
from typing import Union

import numpy as np

from msprime import _msprime

__version__ = "undefined"
try:
    from . import _version

    __version__ = _version.version
except ImportError:
    pass


# Make sure the GSL error handler is turned off so that we can be sure that
# we don't abort on errors. This can be reset by using the function
# _msprime.restore_gsl_error_handler(), which will set the error handler to
# the value that it had before this function was called.
_msprime.unset_gsl_error_handler()


# Some machinery here for generating default random seeds. We need a map
# indexed by process ID here because we cannot use a global variable
# to store the state across multiple processes. Copy-on-write semantics
# for child processes means that they inherit the state of the parent
# process, so if we just keep a global variable without indexing by
# PID, child processes will share the same random generator as the
# parent.

_seed_rng_map: Dict[int, random.Random] = {}


def get_seed_rng() -> Union[random.Random, None]:
    return _seed_rng_map.get(os.getpid(), None)


def clear_seed_rng():
    _seed_rng_map.pop(os.getpid(), None)


def get_random_seed() -> int:
    global _seed_rng_map
    pid = os.getpid()
    if pid not in _seed_rng_map:
        # If we don't provide a seed to Random(), Python will seed either
        # from a system source of randomness (i.e., /dev/urandom) or the
        # current time if this is not available. Thus, our seed rng should
        # be unique, even across different processes.
        _seed_rng_map[pid] = random.Random()
    return _seed_rng_map[pid].randint(1, 2 ** 32 - 1)


def set_seed_rng_seed(seed: int):
    """
    Convenience method to let us make unseeded simulations deterministic
    when generating documentation examples.

    DO NOT USE THIS FUNCTION!!!
    """
    global _seed_rng_map
    pid = os.getpid()
    _seed_rng_map[pid] = random.Random(seed)


def isinteger(value: Any) -> bool:
    """
    Returns True if the specified value can be converted losslessly to an
    integer.
    """
    if isinstance(value, numbers.Number):
        # Mypy doesn't realise we've done an isinstance here.
        return int(value) == float(value)  # type: ignore
    return False


def _parse_flag(value: Any, *, default: bool) -> bool:
    """
    Parses a boolean flag, which can be either True, False, or None.
    If the input value is None, return the default. Otherwise,
    check that the input value is a bool.

    Note that we do *not* cast to a bool as this would accept
    truthy values like the empty list, etc. In this case None
    would be converted to False, potentially conflicting with
    the default value.
    """
    assert isinstance(default, bool)
    if value is None:
        return default
    if not isinstance(value, bool):
        raise TypeError("Boolean flag must be True, False, or None (the default value)")
    return value


@dataclasses.dataclass
class TableEntry:
    data: str
    extra: Union[str, None] = None

    def as_html(self):
        ret = "<td"
        if self.extra is not None:
            wrapped = textwrap.fill(self.extra, 80)
            ret += f" title='{wrapped}'"
        ret += f">{self.data}</td>"
        return ret

    def as_text(self):
        return self.data


def html_table(
    caption: str, column_titles: List[str], data: List[List[Union[TableEntry, str]]]
):
    """
    Returns a HTML table formatted with the specified data.
    """
    header = "".join(f"<th>{col_title}</th>" for col_title in column_titles)
    rows = ""
    for row_data in data:
        assert len(column_titles) == len(row_data)
        row = ""
        for item in row_data:
            if not isinstance(item, TableEntry):
                item = TableEntry(item)
            row += item.as_html()
        rows += f"<tr>{row}</tr>"
    s = (
        "<div>"
        """<style scoped="">
            .tskit-table thead tr th:only-of-type {vertical-align: middle;}
            .tskit-table thead tr th {text-align: center;vertical-align: top;}
            .tskit-table tbody td {text-align: right;padding: 0.5em 0.5em;}
            .tskit-table tbody th {padding: 0.5em 0.5em;}
        </style>"""
        f"<b>{caption}</b>"
        '<table border="1" class="tskit-table">'
        "<thead>"
        "<tr>" + header + "</tr>"
        "</thead>"
        "<tbody>" + rows + "</tbody>"
        "</table>"
        "</div>"
    )
    return s


def _text_table_row(data, alignments, widths):
    num_lines = max(len(item) for item in data)
    for item in data:
        assert isinstance(item, list)
        item.extend([""] * (num_lines - len(item)))
        assert len(item) == num_lines
    s = ""
    for line in range(num_lines):
        out_line = "│"
        for value, align, width in zip(data, alignments, widths):
            out_line += f"{value[line]:{align}{width - 1}}│"
        out_line += "\n"
        s += out_line
    return s


def text_table(
    caption: str,
    column_titles: List[List[str]],
    column_alignments: List[str],
    data: List[List[List[str]]],
    internal_hlines=False,
):
    """
    Returns a text table formatted with the specified data. Column alignments
    should be values used in Python's string formatting mini-language.

    Each item in the table should be a *list* of strings, which are the lines
    of text to be shown in that table cell.
    """
    N = len(column_titles)
    assert len(column_alignments) == N
    widths = np.array([len(title) for title in column_titles], dtype=int)
    for row in data + [column_titles]:
        assert N == len(row)
        for j in range(N):
            widths[j] = max(widths[j], max([len(line) for line in row[j]], default=0))
    widths += 3

    hline = "─" * (sum(widths) - 1)
    internal_hline = "┈" * (sum(widths) - 1)
    out = f"{caption}\n"
    out += f"┌{hline}┐\n"
    out += f"{_text_table_row(column_titles, column_alignments, widths)}"
    out += f"├{hline}┤\n"
    for j, split_row in enumerate(data):
        out += f"{_text_table_row(split_row, column_alignments, widths)}"
        if internal_hlines and j < len(data) - 1:
            out += f"│{internal_hline}│\n"
    out += f"└{hline}┘\n"
    return out
