# -*- coding: utf-8 -*-
# ------------------------------------------------------------------------------
# Description:  Asses the Autocorrelation Coefficient 
# Author:       benedikt.strahm@ilek.uni-stuttgart.de
# Created:      2018-12-06
# Execution:    Import functions / collections (from pyLEK.helpers import util)
#               or executing from command line (py filename.py)
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Sources
# ------------------------------------------------------------------------------
# Sockel, Helmut (1984): Aerodynamik der Bauwerke. Wiesbaden, 
# s.l.: Vieweg+Teubner Verlag., p.21, p. 89
# ------------------------------------------------------------------------------
# Libraries
# ------------------------------------------------------------------------------
import numpy as np

import pyLEK.plotters.plot2D as plt

# ------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------


def calculateCorrelationCoeff(u_prime):
    # Get the size of the array
    sp = u_prime.shape
    nx, ny, nz = sp[0], sp[1], sp[2]

    # Calculate the correlation coefficient
    Ruu = np.zeros((nx, ny, nz))

    for j in range(0, ny):
        for k in range(0, nz):
            u_rever = np.array([u_prime[i, j, k] for i in range(nx-1, -1, -1)])
            u_stack = np.hstack((u_prime[:, j, k], u_prime[:-1, j, k]))
            Ruu[:, j, k] = np.convolve(u_rever, u_stack, mode='valid') / nx

    # For smoothness average over lines in y-direction
    Ruu = np.mean(Ruu, axis=1)

    # Normalization
    Ruu[:, :] = Ruu[:, :] / Ruu[0, :]

    return Ruu


def plotCorrelationCoeff(Ruu, dx, nz):
    # Setting up the x-Axis
    sp = Ruu.shape
    nx = sp[0]
    x = np.linspace(0, (nx-1)*dx, nx)

    # Setting up data to be plotted
    x = [x]
    y = [Ruu[:, nz]]

    # Setting up the plot
    xlabel = "Delta x"
    ylabel = "Ruu (Delta x)"
    title = "Correlation Coefficient"
    legend = ["Ruu at nz = " + str(nz)]

    dir_fileName = "Correlation_Coefficient_Ruu"

    style_dict = {"savefig.format": "pdf", "savefig.format": "pdf"}

    # Plotting
    plt.plot2D(x, y, xlabel=xlabel, ylabel=ylabel, title=title, legend=legend,
               dir_fileName=dir_fileName, style_dict=style_dict, 
               colorScheme='Monochrome', variation='color',
               savePlt=True, showPlt=True,)


def main():
    # --- Input data ---#
    # A three dimensional array containing velocity fluctuations
    # in the order u(x_i, y_i, z_i)
    u_prime =

    # The spacing in time or space of the datapoints
    dx =

    # Choosing at which reference height the
    # correlation coefficient should be plotted
    nz =

    # --- Calculation --#
    # Will result in a form Ruu(x,z)
    Ruu = calculateCorrelationCoeff(u_prime)

    # ---- Plotting ----#
    plotCorrelationCoeff(Ruu, dx, nz)


if __name__ == '__main__':
    main()
