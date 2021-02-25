# -*- coding: utf-8 -*-
# ------------------------------------------------------------------------------
# Description:  Evaluate & plot time series of (global) forces
# Author:       benedikt.strahm@ilek.uni-stuttgart.de
# Created:      2018-01-06
# Execution:    Import functions / collections (from cfdPostProcessing.forces import util)
#               or executing from command line (py filename.py)
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Sources
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Libraries
# ------------------------------------------------------------------------------
import numpy as np
import heapq as hq
import sys
import os
import re

import pyLEK.plotters.plot2D as plt
import pyLEK.helpers.txtEditor as txt

# ------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------


def calcStatistics(Fi):
    # --- Calculation --#
    meanFi = np.mean(Fi)
    stdFi = np.std(Fi)
    rmsFi = np.sqrt(np.mean(np.square(Fi)))
    maxFi = np.max(Fi)
    minFi = np.min(Fi)

    nAvMax = [10, 20, 50, 100]
    maxAvFi = {}
    minAvFi = {}

    for n in nAvMax:
        maxAvFi[n] = np.average(hq.nlargest(n, Fi))
        minAvFi[n] = np.average(hq.nsmallest(n, Fi))

    return meanFi, stdFi, rmsFi, maxFi, minFi, nAvMax, maxAvFi, minAvFi


def writeTimeseries(Fi):
    # --- Calculation --#
    meanFi, stdFi, rmsFi, maxFi, minFi, nAvMax, maxAvFi, minAvFi = calcStatistics(
        Fi)

    comp = "Fx"
    units = "MN"
    dir_fileName = "Timeseries_Baseshear_Fx.txt"

    # Change to current file location
    os.chdir(os.path.dirname(sys.argv[0]))

    # ---- Write to file ----#
    txt.writeToTxt(dir_fileName, "--------------------------")
    txt.writeToTxt(dir_fileName, "Baseforce " + comp + " (" + units + ")")
    txt.writeToTxt(dir_fileName, "--------------------------")
    txt.writeToTxt(dir_fileName, comp + "Mean: " + '{: 10.3f}'.format(meanFi))
    txt.writeToTxt(dir_fileName, comp + "RMS:  " + '{: 10.3f}'.format(rmsFi))
    txt.writeToTxt(dir_fileName, "--------------------------")
    txt.writeToTxt(dir_fileName, comp + "Max:  " + '{: 10.3f}'.format(maxFi))
    txt.writeToTxt(dir_fileName, comp + "Min:  " + '{: 10.3f}'.format(minFi))
    txt.writeToTxt(dir_fileName, comp + "Std:  " + '{: 10.3f}'.format(stdFi))
    txt.writeToTxt(dir_fileName, "--------------------------")

    for n in nAvMax:
        maxAv = maxAvFi[n]
        txt.writeToTxt(dir_fileName, comp + "Max" + str(n) +
                       "Avg:" + '{: 10.3f}'.format(maxAv))
        minAv = minAvFi[n]
        txt.writeToTxt(dir_fileName, comp + "Min" + str(n) +
                       "Avg:" + '{: 10.3f}'.format(minAv))

    txt.writeToTxt(dir_fileName, "--------------------------")


def plotTimeseries(Ti, Fi):
    # --- Calculation --#
    meanFi, stdFi, _rmsFi, _maxFi, _minFi, _nAvMax, _maxAvFi, _minAvFi = calcStatistics(
        Fi)

    hLines = [meanFi, meanFi + stdFi, meanFi - stdFi]
    hTexts = ["\$F_{x,mean}\$", "\$F_{x,mean+std}\$", "\$F_{x,mean-std}\$"]

    # ---- Plotting ----#
    xlabel = "$Time$ $[s]$"
    ylabel = "\$F_{x} [MN]\$"
    title = "Time Series of Base Shear F_{x}"

    # Change to current file location
    os.chdir(os.path.dirname(sys.argv[0]))
    dir_fileName = "Timeseries_Baseshear_Fx"

    xlim = []
    ylim = []

    style_dict = {"lines.linewidth": "0.75", "savefig.format": "svg"}

    plt.plot2D(Ti, Fi, xlabel=xlabel, ylabel=ylabel, title=title,
               dir_fileName=dir_fileName, hLines=hLines, hTexts=hTexts, xlim=xlim, ylim=ylim,
               style_dict=style_dict, colorScheme='UniS', variation='color',
               savePlt=True, showPlt=True)


def importForces():

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

 #1_conv_ref0       || 148-197 || 4900
 #1_conv_ref1       || 168-250 || 8200
 #1_conv_ref2       || 131-250 || 11900
 #2_height8     
 #2_height12        || 0-84    || 8435
 #2_height15
    casefile = "1_conv_ref1"
    pipefile = open(
        '/media/dani/linuxHDD/openfoam/simpleFoam/testing/{}/postProcessing/forces/168/force.dat'.format(casefile), 'r')


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

    t = np.round(t,2)
    fpy = np.array(fpy)

    return fpy


def main():
    # --- Input data ---#
    # List containing time series of forces
    # in the order Fi
    Fi = importForces()  / (10 **6)           # Convert to MN

    # Time stepping
    sT = 168    
    eT = 250
    nT = 8200
    dT = 0.01
    Ti = np.linspace(sT, eT, nT)
 #1_conv_ref0       || 148-197 || 4900
 #1_conv_ref1       || 168-250 || 8200
 #1_conv_ref2       || 131-250 || 11900
 #2_height8     
 #2_height12        || 0-84    || 8435
 #2_height15

    # Specify Cut-Off time
    # sT = int(input("Start Time: "))
    # eT = int(input("End Time:   "))

    # Cut the array to desired range
    # Ti = Ti[int(sT/dT):int(eT/dT)]
    # Fi = Fi[int(sT/dT):int(eT/dT)]

    plotTimeseries(Ti, Fi)
    writeTimeseries(Fi)


if __name__ == '__main__':
    main()
    
