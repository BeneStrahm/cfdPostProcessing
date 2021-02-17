# -*- coding: utf-8 -*-
# ------------------------------------------------------------------------------
# Description:  Asses the Turbulent Kinetic Energy Spectra (PSD)
# Author:       benedikt.strahm@ilek.uni-stuttgart.de
# Created:      2018-01-29
# Execution:    Import functions / collections (from pyLEK.helpers import util)
#               or executing from command line (py filename.py)
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Sources
# ------------------------------------------------------------------------------
# Sockel, Helmut (1984): Aerodynamik der Bauwerke. Wiesbaden,
# s.l.: Vieweg+Teubner Verlag., p.20, p. 97
# ------------------------------------------------------------------------------
# Libraries
# ------------------------------------------------------------------------------
import numpy as np
from scipy.optimize import curve_fit

import pyLEK.plotters.plot2D as plt

# ------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------


def fit_func(n, a, b, beta):

    # Predefinitions:
    alpha = 1
    gamma = 1
    c = 1

    return (a * n ** gamma) / (c + b * n ** alpha) ** beta


def calculateKineticEnergySpecta(u_prime, u_bulk, u_tau, z_ref, u_ref, dx):
    # Get the size of the array
    sp = u_prime.shape
    nx, ny, nz = sp[0], sp[1], sp[2]

    # N is half the length of the output of the FFT (Using symmetrie)
    N = nx//2 + 1

    # Convert to Time domain, see Taylors frozen turbulence hyptothesis
    dT = dx / u_bulk

    # Nyquist frequency / wave number
    fNyq = 1 / (2 * dT)                                  # Nyquist frequency
    # kNyq    = 1 / (2 * dx)                                  # Nyquist wave number

    # For explanation see https://www.cbcity.de/die-fft-mit-python-einfach-erklaert
    # Determine frequencies resulting from FFT

    # Space domain:
    # -------------
    # kx      = np.linspace(0, kNyq, N, endpoint=True)            # Wave numbers
    # mu      = kx * 2.0 * np.pi                                  # Circular wave number
    # kx_red   = kx * z_ref                                        # Reduced Wave Number

    # Time domain
    # -------------
    fx = np.linspace(0, fNyq, N, endpoint=True)            # Frequency
    omega = fx * 2.0 * np.pi                                  # Circular frequency
    fx_red = fx * z_ref / u_ref                                # Reduced frequency

    # Empty power spectral density
    # Due to symmetrie, only the first half of the FFT is of interest
    PSD_u_prime = np.zeros(N)

    # Loop in z-direction
    for k in range(0, nz):
        # Loop in y-direction
        for j in range(0, ny):

            # Calculate the wind spectrum
            PSD_u_prime_ = abs(np.fft.fft(u_prime[:, j, k])[:N])

            # Get the Power Spectral Density
            PSD_u_prime_ = PSD_u_prime_**2

            # Scaling Factor
            PSD_u_prime_ = PSD_u_prime_ / nx

            # Scaling Factor
            # PSD_u_prime_ = PSD_u_prime_ * dx                    # Space domain
            PSD_u_prime_ = PSD_u_prime_ * dT                    # Time domain

            # Normalize by friction velocity and frequency
            # PSD_u_prime_ = PSD_u_prime_ * kx / (utau ** 2)      # Space domain
            PSD_u_prime_ = PSD_u_prime_ * fx / (u_tau ** 2)     # Time domain

            # Add to average
            PSD_u_prime += 1 / (ny * nz) * PSD_u_prime_

    return PSD_u_prime, fx_red

    # END
    # ----------------------------------------------------------------------


def plot(PSD_u_prime_measured, fx_red):
    # --- To compare with exact spectra --- #
    # Exact Kaimal Spectra, see "Kaimal, Wyngaard, Izumi, Cote - Spectral
    # characteristics of surface layer turbulence", Quarterly Journal of the Royal
    # Meteorological Society 98 (1972) or "Cheynet, Jakobsen, Obhrai - Spectral characteristics
    # of surface-layer turbulence in the North Sea" (2017)

    def PSD_u_prime_exact(n): return (a * n ** gamma) / \
        (c + b * n ** alpha) ** beta

    # Kaimal Spectra:
    a = 1 * 52.5                      # Factor 1 or 2 for one- or two-sided specta
    b = 33
    alpha = 1
    beta = 5/3
    gamma = 1
    c = 1

    # --- Fit function to measured spectra --- #
    # Fit a function to the values
    popt, pcov = curve_fit(fit_func, fx_red, PSD_u_prime_measured)
    PSD_u_prime_fit = fit_func(fx_red, *popt)

    # Setting up data to be plotted
    x = [fx_red]
    y = [PSD_u_prime_exact(fx_red), PSD_u_prime_measured, PSD_u_prime_fit]

    # Setting up the plot
    xlabel = "fx,red"
    ylabel = "f x PSD / (u_tau ** 2)"
    title = "Turbulent kinetic energy spectrum"
    legend = ["exact", "measured", "fit"]

    dir_fileName = "Turbulenc_kinetic_energy_spectrum_Su"

    style_dict = {"savefig.format": "pdf", "savefig.format": "pdf"}

    # Plotting
    plt.plot2D(x, y, xlabel=xlabel, ylabel=ylabel, title=title, legend=legend,
               xscale='log', yscale='log',
               dir_fileName=dir_fileName, style_dict=style_dict,
               colorScheme='Monochrome', variation='color',
               savePlt=True, showPlt=True)


def main():
    # --- Input data ---#
    # A three dimensional array containing velocity fluctuations
    # in the order u(x_i, y_i, z_i)
    u_prime =

    # The spacing of the datapoints
    dx =

    # Bulk velocity
    u_bulk =

    # Friction velocity
    u_tau =

    # Characteristic height and velocity for reduced frequency
    u_ref =
    z_ref =

    # --- Calculation --#
    # Will result in a form Iuu(z)
    Su, fx_red = calculateKineticEnergySpecta(
        u_prime, u_bulk, u_tau, z_ref, u_ref, dx)

    # ---- Plotting ----#
    plot(Su, fx_red)


if __name__ == '__main__':
    main()
