# -*- coding: utf-8 -*-
# ------------------------------------------------------------------------------
# Description:  ** Add here short description **
# Author:       ** Add here author's e-mail adress **
# Created:      ** Add here the date of creation **
# Execution:    Import functions / collections (from helpers import util)
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Sources
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Libraries
# ------------------------------------------------------------------------------

import re
import numpy as np

# ------------------------------------------------------------------------------
# Functions


def importForces(fname):

    forceRegex = r"([0-9.Ee\-+]+)\s+\(+([0-9.Ee\-+]+)\s([0-9.Ee\-+]+)\s([0-9.Ee\-+]+)+\)+\s+\(+([0-9.Ee\-+]+)\s([0-9.Ee\-+]+)\s([0-9.Ee\-+]+)+\)+\s+\(+([0-9.Ee\-+]+)\s([0-9.Ee\-+]+)\s([0-9.Ee\-+]+)+\)"

    t = []
    ftotx = []
    ftoty = []
    ftotz = []  # Porous
    fpx = []
    fpy = []
    fpz = []  # Pressure
    fvx = []
    fvy = []
    fvz = []  # Viscous

    pipefile = open(fname, 'r')

    lines = pipefile.readlines()

    for line in lines:
        match = re.search(forceRegex, line)
        if match:
            t.append(float(match.group(1)))
            ftotx.append(float(match.group(2)))
            ftoty.append(float(match.group(3)))
            ftotz.append(float(match.group(4)))
            fpx.append(float(match.group(5)))
            fpy.append(float(match.group(6)))
            fpz.append(float(match.group(7)))
            fvx.append(float(match.group(8)))
            fvy.append(float(match.group(9)))
            fvz.append(float(match.group(10)))

    # t = np.round(t, 2)
    fpy = np.array(fpy)

    return fpy


def sample():                                       # Sample of the function
    pass


if __name__ == "__main__":
    sample()
