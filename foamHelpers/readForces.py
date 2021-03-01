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


# interpoliert forces in x und y richtung für delta 0.01 
    tmin = min(t)
    tmax = max(t)
    t_interp = np.round(np.linspace(tmin,tmax,tmax/0.01,endpoint=True),2)
    fpylen = len(fpy)
    fpy_interp = np.interp(t_interp, t, fpy)
    fpx_interp = np.interp(t_interp, t, fpx)


    fpy = np.array(fpy_interp)
    return t_interp, fpy_interp, fpx_interp





if __name__ == "__main__":
    importForces()
