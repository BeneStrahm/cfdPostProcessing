# ------------------------------------------------------------------------------
# Project:      Wind tunnel postprocessing
# Description:  Main file
# Author:       benedikt.strahm@ilek.uni-stuttgart.de
# Created:      2020-09-15
# Execution:    Drag & drop file on script
#               https://mindlesstechnology.wordpress.com/2008/03/29/make-python-scripts-droppable-in-windows/
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Libraries
# ------------------------------------------------------------------------------  

import sys           
import numpy as np      
import heapq as hq       
import os                                

# ------------------------------------------------------------------------------
# Imported functions
# ------------------------------------------------------------------------------

import cfdPostProcessing.WindTunnelPostprocessing.aeroForces as aeroForces
import cfdPostProcessing.WindTunnelPostprocessing.modelProp as modelProp
import cfdPostProcessing.WindTunnelPostprocessing.scaling as scaling
import pyLEK.helpers.txtEditor as txt
import pyLEK.plotters.plot2D as plt

from pyLEK.helpers.pyExtras import getKeyList
from pyLEK.helpers.filemanager import delFilesInFolder

# ------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------     

def calcStatistics(BF_p):
    # --- Calculation --#
    meanBF = np.mean(BF_p)
    stdBF = np.std(BF_p)
    rmsBF = np.sqrt(np.mean(np.square(BF_p)))
    maxBF = np.max(BF_p)
    minBF = np.min(BF_p)

    nAvMax = [10, 20, 50, 100]
    maxAvBF = {}
    minAvBF = {}

    for n in nAvMax:
        maxAvBF[n] = np.average(hq.nlargest(n, BF_p))
        minAvBF[n] = np.average(hq.nsmallest(n, BF_p))

    return meanBF, stdBF, rmsBF, maxBF, minBF, nAvMax, maxAvBF, minAvBF

def writeTimeseries(BF_p):
    # --- Calculation --#
    meanBF, stdBF, rmsBF, maxBF, minBF, nAvMax, maxAvBF, minAvBF = calcStatistics(
        BF_p)

    units = "MN"
    dir_fileName = "WindTunnel.txt"

    # Change to current file location
    os.chdir(os.path.dirname(sys.argv[0]))

    # ---- Write to file ----#
    txt.writeToTxt(dir_fileName, "--------------------------")
    txt.writeToTxt(dir_fileName, "Baseforce "" (" + units + ")")
    txt.writeToTxt(dir_fileName, "--------------------------")
    txt.writeToTxt(dir_fileName, "Mean: " + '{: 10.3f}'.format(meanBF))
    txt.writeToTxt(dir_fileName, "RMS:  " + '{: 10.3f}'.format(rmsBF))
    txt.writeToTxt(dir_fileName, "--------------------------")
    txt.writeToTxt(dir_fileName, "Max:  " + '{: 10.3f}'.format(maxBF))
    txt.writeToTxt(dir_fileName, "Min:  " + '{: 10.3f}'.format(minBF))
    txt.writeToTxt(dir_fileName, "Std:  " + '{: 10.3f}'.format(stdBF))
    txt.writeToTxt(dir_fileName, "--------------------------")

    for n in nAvMax:
        maxAv = maxAvBF[n]
        txt.writeToTxt(dir_fileName, "Max" + str(n) +
                       "Avg:" + '{: 10.3f}'.format(maxAv))
        minAv = minAvBF[n]
        txt.writeToTxt(dir_fileName, "Min" + str(n) +
                       "Avg:" + '{: 10.3f}'.format(minAv))

    txt.writeToTxt(dir_fileName, "--------------------------")


def main():
    # Get name of input file
    # fname = sys.argv [1]
    fname = "/media/dani/linuxHDD/visualStudio/cfdPostProcessing/WindTunnelPostprocessing/T115_6/T115_6_000.mat"

    # Clean up results folder
    # delFilesInFolder('T115_6/results')

    # Full scale building properties
    uH_f    = 36             # m/s       // Wind speed at z = H (50yr)
    H_f     = 160               # m         // Building height
    B       = 32                # m         // Building width
    dns     = ['D', 'L']        #           // Directions (Drag/Lifts)
 

    # @DL: Folgende Variablen ignorieren, 
    # Diese wären nur wichtig für Strukturanalyse
    nF      = 32                #           // Number of floors
    nM      = 4                 #           // Number of modules
    b       = 16                # m         // Core wall thickness    
    D       = 0.02              # %         // Damping
    I       = 477.924           # m4        // Starting value
    E       = 28900 * 10 ** 3   # kN/m2     // E-Modulus
    mue     = 30473 / H_f       # t/m       // Mass distribution

    # Iterate over both directions (D/L)
    for dn in dns:    
        # Load wind tunnel model properties, TPU Database files
        wtModelProp = modelProp.wtModelProp(fname)
        wtModelProp.loadTPUModelProp()

        # Initialize building model properties
        buildProp = modelProp.buildProp(H_f, B, nF, nM, dn, E, I, D, uH_f)
        
        # Load aerodynamic forces in model scale
        wtModelAeroForces = aeroForces.wtModelAeroForces()
        wtModelAeroForces.loadTPUModelForces(wtModelProp, buildProp)

        # Calculate scaling factors
        scalingFactors = scaling.scalingFactors(wtModelProp, buildProp)

        # Scale wind tunnel model
        buildProp.scaleBuildProp(wtModelProp, scalingFactors)

        # Scale forces
        buildAeroForces = aeroForces.buildAeroForces(scalingFactors, wtModelAeroForces)

        # All forces are in kN
        print("Direction: " + dn)
        print("Mean Base Force: " + '{:02.3f}'.format(np.mean(buildAeroForces.BF_p)))
        print("Min Base Force: " + '{:02.3f}'.format(np.min(buildAeroForces.BF_p)))
        print("Max Base Force: " + '{:02.3f}'.format(np.max(buildAeroForces.BF_p)))
        print("Std Base Force: " + '{:02.3f}'.format(np.std(buildAeroForces.BF_p)))
        
        meanBF, stdBF, rmsBF, maxBF, minBF, nAvMax, maxAvBF, minAvBF = calcStatistics(buildAeroForces.BF_p)
        hLines = [meanBF, meanBF + stdBF, meanBF - stdBF]
        hTexts = ["\$F_{x,mean}\$", "\$F_{x,mean+std}\$", "\$F_{x,mean-std}\$"]




        # Time series of base forces
        t = np.linspace(0, buildProp.nT*buildProp.dT, buildProp.nT)
        F = buildAeroForces.BF_p
        writeTimeseries(F)
        
        style_dict = {"lines.linewidth": "0.5", "savefig.format": "svg"}
        plt.plot2D(t, F, hLines=hLines,   style_dict=style_dict, xlabel="Time [s]", ylabel= "Base Force " + dn + " [kN]", title="Time series of base forces", showPlt=True)            
        print(stdBF)
            

if __name__ == '__main__':
    main()
    # wait = input("Press Enter to continue.")

