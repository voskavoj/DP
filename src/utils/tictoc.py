"""
    Timers for measuring time of execution of code
"""

import time


class TicToc:
    def __init__(self):
        self.act_time = time.monotonic()

    def tic(self):
        self.act_time = time.monotonic()

    def toc(self, note=""):
        print(note, time.monotonic() - self.act_time)
        self.tic()


TICTOC = TicToc()


def tic():
    global TICTOC
    TICTOC.tic()


def toc(note=""):
    global TICTOC
    TICTOC.toc(note)
