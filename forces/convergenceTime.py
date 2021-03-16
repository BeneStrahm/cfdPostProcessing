# -*- coding: utf-8 -*-
# ------------------------------------------------------------------------------
# Description:  Evaluate & plot convergence time of (global) forces
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

import pyLEK.plotters.plot2D as plt

from foamHelpers import readForces

# ------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------
1345

def statisticsOverTimeSeries(data):
    meanBF = []
    rmsBF  = []
    stdBF  = []

    for i in range(0, len(data)):
        # Get current mean
        meanBF.append(np.mean(data[0:i+1]))     

        # Get current rms
        rmsBF.append(np.sqrt(np.mean(np.square(data[0:i+1])))) 

        # Get current mean
        stdBF.append(np.std(data[0:i+1]))

    return meanBF, rmsBF, stdBF

def main(dObject):
    # Get ouput directory
    outDir = '/media/dani/linuxHDD/openfoam/simpleFoam/testing/'
    fname0 = '/media/dani/linuxHDD/openfoam/simpleFoam/testing/1_conv_ref0/postProcessing/forces/0/force.dat'
    interpForces = readForces.importForces(fname0)
    

        # Get data to plot
    sT = min(interpForces[0])   
    eT = max(interpForces[0])
    dT = interpForces[0][1]-interpForces[0][0]
    nT = len(interpForces[0]*dT)
        
    T = np.linspace(sT, eT, nT)

    BF = interpForces[2]/ (10 ** 6) 

    # Specify Cut-Off time
    print("Specify time range to plot")
    sT = int(input("Start Time: "))
    eT = int(input("End Time:   "))
    print("Specify time to start evaluation")
    sT_S = int(input("Start Time: "))

    for comp in BF:  
        
        # Cut the array to desired range 
        x = T[int(sT_S/dT):int(eT/dT)]
        y = BF[comp][int(sT_S/dT):int(eT/dT)]

        # Filter anomalies if desired
        # flagFilter = input("Shall data for " + comp + " be filtered? (Y/N): ")
    
        # if flagFilter == "Y":
        #     y = filterDataWithStd(y, comp, dT, sT)

        # Calculate statistic properties     
        meanBF, rmsBF, stdBF  = statisticsOverTimeSeries(y)

        x = [x]

        # Plot the mean / time  
        y = [meanBF]

        xlabel = r"$Time [s]$"
        if "F" in comp:
            ylabel = r"$Mean F_{} [MN]$".format(comp[1])
            legend = [r"$Mean F_{}$".format(comp[1])]
            title = "Mean Base Shear"
            app = comp + "_" + str(int(round(np.min(x)))) + "_" + str(int(round(np.max(x))))
            dir_fileName = outDir + "Mean_BaseShear_" + app

        elif "M" in comp:
            ylabel = r"$Mean M_{} [MNm]$".format(comp[1])
            legend = [r"$Mean M_{}$".format(comp[1])]
            title = "Mean Base Moment"
            app = comp + "_" + str(int(round(np.min(x)))) + "_" + str(int(round(np.max(x))))
            dir_fileName = outDir + "Mean_BaseMoment_" + app

        xlim = []
        ylim = []

        style_dict = {"lines.linewidth":1.0, "figure.figsize" : "11.0, 4.8"}

        plt.plot2D(x, y, xlabel, ylabel, title, legend, dir_fileName, 
                    xlim=xlim, ylim=ylim, xscale='linear', yscale='linear',
                    style_dict=style_dict, colorScheme='TUM', variation='color')

        # Plot the root-mean-square / time 
        y = [rmsBF]

        xlabel = r"$Time [s]$"
        if "F" in comp:
            ylabel = r"$RMS F_{} [MN]$".format(comp[1])
            legend = [r"$RMS F_{}$".format(comp[1])]
            title = "Root-Mean-Square Base Shear"
            app = comp + "_" + str(int(round(np.min(x)))) + "_" + str(int(round(np.max(x))))
            dir_fileName = outDir + "RMS_BaseShear_" + app

        elif "M" in comp:
            ylabel = r"$RMS M_{} [MNm]$".format(comp[1])
            legend = [r"$RMS M_{}$".format(comp[1])]
            title = "Root-Mean-Square Base Moment"
            app = comp + "_" + str(int(round(np.min(x)))) + "_" + str(int(round(np.max(x))))
            dir_fileName = outDir + "RMS_BaseMoment_" + app

        xlim = []
        ylim = []

        style_dict = {"lines.linewidth":1.0, "figure.figsize" : "11.0, 4.8"}

        plt.plot2D(x, y, xlabel, ylabel, title, legend, dir_fileName, 
                    xlim=xlim, ylim=ylim, xscale='linear', yscale='linear',
                    style_dict=style_dict, colorScheme='TUM', variation='color')

        # Plot the standart deviation / time 
        y = [stdBF]

        xlabel = r"$Time [s]$"
        if "F" in comp:
            ylabel = r"$STD F_{} [MN]$".format(comp[1])
            legend = [r"$STD F_{}$".format(comp[1])]
            title = "Standard Deviation Base Shear"
            app = comp + "_" + str(int(round(np.min(x)))) + "_" + str(int(round(np.max(x))))
            dir_fileName = outDir + "STD_BaseShear_" + app

        elif "M" in comp:
            ylabel = r"$STD M_{} [MNm]$".format(comp[1])
            legend = [r"$STD M_{}$".format(comp[1])]
            title = "Standard Deviation Base Moment"
            app = comp + "_" + str(int(round(np.min(x)))) + "_" + str(int(round(np.max(x))))
            dir_fileName = outDir + "STD_BaseMoment_" + app

        xlim = []
        ylim = []

        style_dict = {"lines.linewidth":1.0, "figure.figsize" : "11.0, 4.8"}

        plt.plot2D(x, y, xlabel, ylabel, title, legend, dir_fileName, 
                    xlim=xlim, ylim=ylim, xscale='linear', yscale='linear',
                    style_dict=style_dict, colorScheme='TUM', variation='color')

    print("\n")



def main():
    # --- Input data ---#
    fname = '/media/dani/linuxHDD/openfoam/simpleFoam/testing/17_v1Fine/postProcessing/forces/0/force.dat'
    Fi = readForces.importForces(fname) / (10 ** 6)              # Convert to MN    
    # --- Calculation --#

    # ---- Plotting ----#


if __name__ == '__main__':
    main()
