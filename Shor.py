# Shor's Algorithm - Introductory Quantum Algorithm
# Start on August 6th 2018
# Implemented by Jay Whaley

from __future__ import print_function
from halo import Halo
import cirq
import random
import sys
import numpy as np


# Standard Error Stream Printing: from https://stackoverflow.com/a/14981125/7476183
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)
    pass


pass


def compute_gcd(x, y):
    while y:
        x, y = y, x % y

        return x
    pass


def classical_part(n):
    a = random.randint()
    gcd = compute_gcd(a, n)
    if gcd != 1:
        return
    else:

        pass
    pass


def period_finding_subroutine(n):
    grid_length = 4
    pass


pass
