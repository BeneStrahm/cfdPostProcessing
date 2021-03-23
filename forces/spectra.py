# -*- coding: utf-8 -*-
# ------------------------------------------------------------------------------
# Description:  Evaluate the (global) force specta
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
from matplotlib.pyplot import yscale
import numpy as np
import sys
import os
import re
import pyLEK.plotters.plot2D as plt

from cfdPostProcessing.foamHelpers import readForces

# ------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------


def calculateSpectra(y, dT):
    # Get shape of truncated forces
    sp = y.shape
    nT = sp[0]

    # N is half the length of the output of the FFT (Using symmetrie)
    N = nT//2 + 1            # // -> int

    # Calculate the Nyquist frequency
    fNyq = 1 / (2 * dT)                                  # Nyquist frequency

    # Empty power spectral density
    # Due to symmetrie, only the first half of the FFT is of interest
    Sa = np.zeros(N)

    # For explanation see https://www.cbcity.de/die-fft-mit-python-einfach-erklaert
    # Determine frequencies resulting from FFT
    # https://github.com/pyNFFT/pyNFFT for non uniform samples
    # Time domain
    # -------------
    f = abs(np.fft.fftfreq(nT, dT)[:N])                         # Frequency
    # f      = np.linspace(0, fNyq, N, endpoint=True)          # Same as above

    # Calculate the force spectrum
    Sa = abs(np.fft.fft(y)[:N])

    # Get the Power Spectral Density
    Sa = Sa**2

    # Scaling Factor
    Sa = Sa / nT

    # Scaling Factor
    Sa = Sa * dT

    # Normalize by standart deviation and frequency
    Sa = Sa * f / (np.std(y) ** 2)

    return f, Sa


def plotSpectra(f, Sa):         #comp
    # Function fit
    # # Set up logspace
    # fx_fit = np.logspace(np.log10(fx[1]), np.log10(fx[-1]), num=1000)

    # # Interpolate to grid
    # Sa_fit = griddata(fx, Sa, fx_fit, method='cubic')

    # Fit a function to the values
    # Sa_fit = savgol_filter(Sa, 51, 3)      # window size 51, polynomial order 3

    # Setting up the plot
    x = [f]
    y = [Sa]

    #maximum peak frequency




    max_fq = f[np.argmax(Sa)]

    # ---- Plotting ----#
    xlabel = "f"
    ylabel = "PSD " #(" + comp + ")
    title = "Power Spectral Density"
    vLines = [max_fq]
    vTexts = [r"\$f_{peak}" + '{:03.4f}'.format(max_fq) +r"\$"]
    legend = ["measured"]

    # Change to current file location
    os.chdir(os.path.dirname(sys.argv[0]))
    dir_fileName = "PowerSpectralDensity__Baseshear_Fx"

    xlim = []
    ylim = [10**-3, 10**1.5]

    style_dict = {"savefig.format": "svg"}

    plt.plot2D(x, y, xlabel=xlabel, ylabel=ylabel, title=title, legend=legend, vLines=vLines, vTexts=vTexts,
               dir_fileName=dir_fileName, xlim=xlim, ylim=ylim,
               xscale='log', yscale='log',
               style_dict=style_dict, colorScheme='UniS', variation='color',
               savePlt=True, showPlt=True, savePkl=True,)

def main():
    # --- Input data ---#
    fname = '/media/dani/linuxHDD/openfoam/simpleFoam/testing/1_conv_ref2/postProcessing/forces/0/force.dat'
    interpForces = readForces.importForces(fname)
    y = interpForces[1]             # Convert to MN
    dT = interpForces[0][1]-interpForces[0][0]
        # Specify Cut-Off time
    # sT = int(input("Start Time: 100 "))
    # eT = int(input("End Time: 150   "))
    sT = 100
    eT = 200


    # Cut the array to desired range
    y = y[int(sT/dT):int(eT/dT)]
    f, Sa = calculateSpectra(y,dT)
    plotSpectra(f, Sa)
    print(f[np.argmax(Sa)])


if __name__ == '__main__':
    main()
