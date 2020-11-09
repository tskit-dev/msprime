"""
Code to stress the low-level API as much as possible to expose
any memory leaks or error handling issues.
"""
import curses
import os
import random
import resource
import sys
import time
import tracemalloc
from contextlib import redirect_stdout

import pytest


def main(stdscr):
    if len(sys.argv) > 1:
        args = sys.argv[1:]
    else:
        args = ["-n0", "tests/test_lowlevel.py"]

    class StressPlugin:
        def __init__(self):
            self.max_rss = 0
            self.max_rss_iter = 0
            self.min_rss = 1e100
            self.iteration = 0
            self.last_print = time.time()
            self.memory_start = None

        def pytest_sessionstart(self):
            if self.memory_start is None:
                tracemalloc.start()
                self.memory_start = tracemalloc.take_snapshot()

        def pytest_sessionfinish(self):
            memory_current = tracemalloc.take_snapshot()
            rusage = resource.getrusage(resource.RUSAGE_SELF)
            if self.max_rss < rusage.ru_maxrss:
                self.max_rss = rusage.ru_maxrss
                self.max_rss_iter = self.iteration
            if self.min_rss > rusage.ru_maxrss:
                self.min_rss = rusage.ru_maxrss

            # We don't want to flood stdout, so we rate-limit to 1 per second.
            if time.time() - self.last_print > 1:
                stdscr.clear()
                rows, cols = stdscr.getmaxyx()
                stdscr.addstr(
                    0,
                    0,
                    "iter\tRSS\tmin\tmax\tmax@iter"[: cols - 1],
                )
                stdscr.addstr(
                    1,
                    0,
                    "\t".join(
                        map(
                            str,
                            [
                                self.iteration,
                                rusage.ru_maxrss,
                                self.min_rss,
                                self.max_rss,
                                self.max_rss_iter,
                            ],
                        )
                    )[: cols - 1],
                )
                stats = memory_current.compare_to(self.memory_start, "traceback")
                for i, stat in enumerate(stats[: rows - 3], 1):
                    stdscr.addstr(i + 2, 0, str(stat)[: cols - 1])
                self.last_print = time.time()
                stdscr.refresh()
                self.iteration += 1

    plugin = StressPlugin()
    while True:
        # We don't want any random variation in the amount of memory
        # used from test-to-test.
        random.seed(1)
        with open(os.devnull, "w") as devnull:
            with redirect_stdout(devnull):
                result = pytest.main(args, plugins=[plugin])
        if result != 0:
            exit("TESTS FAILED")


if __name__ == "__main__":
    stdscr = curses.initscr()
    curses.noecho()
    curses.cbreak()

    try:
        main(stdscr)
    finally:
        curses.echo()
        curses.nocbreak()
        curses.endwin()
