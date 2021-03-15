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


def importMoments(mname):

    forceRegex = r"([0-9.Ee\-+]+)\s+\(+([0-9.Ee\-+]+)\s([0-9.Ee\-+]+)\s([0-9.Ee\-+]+)+\)+\s+\(+([0-9.Ee\-+]+)\s([0-9.Ee\-+]+)\s([0-9.Ee\-+]+)+\)+\s+\(+([0-9.Ee\-+]+)\s([0-9.Ee\-+]+)\s([0-9.Ee\-+]+)+\)"

    t = []
    mtotx = []
    mtoty = []
    mtotz = []  # Porous
    mpx = []
    mpy = []
    mpz = []  # Pressure
    mvx = []
    mvy = []
    mvz = []  # Viscous

    pipefile = open(mname, 'r')

    lines = pipefile.readlines()

    for line in lines:
        match = re.search(forceRegex, line)
        if match:
            t.append(float(match.group(1)))
            mtotx.append(float(match.group(2)))
            mtoty.append(float(match.group(3)))
            mtotz.append(float(match.group(4)))
            mpx.append(float(match.group(5)))
            mpy.append(float(match.group(6)))
            mpz.append(float(match.group(7)))
            mvx.append(float(match.group(8)))
            mvy.append(float(match.group(9)))
            mvz.append(float(match.group(10)))


# interpoliert forces in x und y richtung für delta 0.01 
    tmin = min(t)
    tmax = max(t)

    t_interp = np.round(np.linspace(tmin,tmax,int(tmax/0.01),endpoint=True),2)
    My = np.interp(t_interp, t, mpy)
    Mx = np.interp(t_interp, t, mpx)
    Mz = np.interp(t_interp, t, mpz)


    return t_interp, My, Mx, Mz





if __name__ == "__main__":
    importMoments()
