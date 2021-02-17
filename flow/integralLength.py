# -*- coding: utf-8 -*-
# ------------------------------------------------------------------------------
# Description:  Asses the Integral Length Scale
# Author:       benedikt.strahm@ilek.uni-stuttgart.de
# Created:      2018-12-06
# Execution:    Import functions / collections (from cfdPostProcessing.flow import util)
#               or executing from command line (py filename.py)
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Sources
# ------------------------------------------------------------------------------
# Sockel, Helmut (1984): Aerodynamik der Bauwerke. Wiesbaden, 
# s.l.: Vieweg+Teubner Verlag., p.22, p. 89
# ------------------------------------------------------------------------------
# Libraries
# ------------------------------------------------------------------------------
import numpy as np
import os, sys

import pyLEK.plotters.plot2D as plt

# ------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------


def calculateIntegralLength(Ruu, dx):   
    # Get the size of the array
    sp = Ruu.shape
    nx, nz = sp[0], sp[1]  

    # Calculate the integral length
    Lux = np.zeros((nx, nz))

    for k in range(0, nz):
        
        # Calculating integral length using numerical integration
        for i in range(0,nx-1):                
            Lux[i+1,k] = Lux[i,k] + (Ruu[i,k]+Ruu[i+1,k]) / 2 * dx

    # Determine maximal values over the height
    Lux_max = np.zeros(nz)

    # As the correlation coefficient can have values < 0 and is symmetric
    # the maximal integral length can be determined by dLux/dx < 0
    for k in range(0,nz):
        i = 0
        while Lux[i,k] < Lux[i+1,k] and i < len(Lux[:,k]) - 2:
            i += 1 
        Lux_max[k] = Lux[i,k]    

    Lux = Lux_max

    return Lux


def plotIntegralLength(Lux, dz):
    # Setting up the y-Axis
    sp  = Lux.shape
    nz  = sp[0]
    z   = np.linspace(0, (nz-1)*dz, nz)

    # Setting up data to be plotted
    x = [Lux]
    y = [z]

    # Setting up the plot
    xlabel  = "Lux(z)"
    ylabel  = "z"
    title   = "Integral length scale"
    legend  = ["Lux(z)"]

    # Change to current file location
    os.chdir(os.path.dirname(sys.argv[0]))
    dir_fileName = "Integral_Length_Luu"

    style_dict = {"savefig.format": "svg"}

    # Plotting
    plt.plot2D(x, y, xlabel=xlabel, ylabel=ylabel, title=title, legend=legend,
               dir_fileName=dir_fileName, style_dict=style_dict, 
               colorScheme='UniS', variation='color',
               savePlt=True, showPlt=True)

def main():
    # --- Input data ---#
    # The autocorrelation coefficient in the form Ruu(x_i, z_i)
    Ruu     = 

    # The spacing in time or space of the datapoints
    dx      = 
    dz      = 

    # --- Calculation --# 
    # Will result in a form Lux(z)
    Lux = calculateIntegralLength(Ruu, dx)

    # ---- Plotting ----#
    plotIntegralLength(Lux, dz)

if __name__ == '__main__':
    main()