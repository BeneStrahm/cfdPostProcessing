# -*- coding: utf-8 -*-
# ------------------------------------------------------------------------------
# Description:  Asses the Turbulence Intensity
# Author:       benedikt.strahm@ilek.uni-stuttgart.de
# Created:      2018-12-06
# Execution:    Import functions / collections (from pyLEK.helpers import util)
#               or executing from command line (py filename.py)
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Sources
# ------------------------------------------------------------------------------
# Sockel, Helmut (1984): Aerodynamik der Bauwerke. Wiesbaden,
# s.l.: Vieweg+Teubner Verlag., p.19, p. 88
# ------------------------------------------------------------------------------
# Libraries
# ------------------------------------------------------------------------------
import numpy as np

import pyLEK.plotters.plot2D as plt

# ------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------


def calculateTurbulenceIntensity(u_prime, dz, kappa, z0, utau):
    # Get the size of the array
    sp = u_prime.shape
    ny, nz = sp[1], sp[2]

    # Get the standart deviation and...
    std_u_prime_ = np.zeros((ny, nz))
    for k in range(0, ny):
        std_u_prime_[k, :] += np.std(u_prime[:, k, :], axis=(0))

    # ...average over the y-direction
    std_u_prime = np.mean(std_u_prime_, axis=0)

    # Create logarhytmic wind profile
    z = np.linspace(0, (nz-1)*dz, nz)
    ulog = np.where(z > 0, utau / kappa * np.log(z/z0), 0)

    # Create the turbulence intensity profiles
    Iuu = std_u_prime / ulog

    return Iuu


def plotTurbulenceIntensity(Iuu, dz):
    # Setting up the y-Axis
    sp = Iuu.shape
    nz = sp[0]
    z = np.linspace(0, (nz-1)*dz, nz)

    # Setting up data to be plotted
    x = [Iuu*100]
    y = [z]

    # Setting up the plot
    xlabel = "Iuu(z)"
    ylabel = "z"
    title = "Turbulence intensity"
    legend = ["Iuu(z)"]

    dir_fileName = "Turbulence_Intensity_Iuu"

    style_dict = {"savefig.format": "pdf", "savefig.format": "pdf"}

    # Plotting
    plt.plot2D(x, y, xlabel=xlabel, ylabel=ylabel, title=title, legend=legend,
               dir_fileName=dir_fileName, style_dict=style_dict,
               colorScheme='Monochrome', variation='color',
               savePlt=True, showPlt=True)


def main():
    # --- Input data ---#
    # A three dimensional array containing velocity fluctuations
    # in the order u(x_i, y_i, z_i)
    u_prime =

    # The spacing of the datapoints
    dz =

    # Parameters for the target wind profile
    kappa =
    z0 =
    utau =

    # --- Calculation --#
    # Will result in a form Iuu(z)
    Iuu = calculateTurbulenceIntensity(u_prime, dz, kappa, z0, utau)

    # ---- Plotting ----#
    plotTurbulenceIntensity(Iuu, dz)


if __name__ == '__main__':
    main()
