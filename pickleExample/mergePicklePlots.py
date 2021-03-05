#!/usr/bin/env python3
# encoding: utf-8

# ------------------------------------------------------------------------------
# Project       : Master Thesis
# Description   : Merge Plots which were saved as .pickle
# Command       : -
# Author        : benedikt.strahm@tum.de
# Created       : 2019-03-06
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Imported functions
# ------------------------------------------------------------------------------

import pickle as pkl
import pyLEK.plotters.plot2D as plot2D
import os,sys
import matplotlib.pyplot as plt

def getData(ax):
    line = ax.lines[0]
    x = line.get_xdata()
    y = line.get_ydata()

    return x, y

def main():

    # Specify the file name
    # fname = str(input("Specify File Name: "))
    fname1 = "Timeseries_Baseshear_Fx_1_conv_ref0"
    fname2 = "Timeseries_Baseshear_Fx_1_conv_ref1"
    fname3 = "Timeseries_Baseshear_Fx_1_conv_ref2"

    # Change to current     
    # Change to current file location
    os.chdir(os.path.dirname(sys.argv[0]))

    # Load the ax-object
    ax1 = pkl.load(open(fname1 + '.pickle', 'rb'))
    ax2 = pkl.load(open(fname2 + '.pickle', 'rb'))
    ax3 = pkl.load(open(fname3 + '.pickle', 'rb'))

    # Get data from the ax-object
    x1, y1 = getData(ax1)
    x2, y2 = getData(ax2)
    x3, y3 = getData(ax3)

    # Clean up
    plt.close('all')

    # Plotting 
    x = [x1+2, x2+1, x3]
    y = [y1, y2, y3]


    plot2D.plot2D(x, y, showPlt=True)

if __name__ == '__main__':
    main()
    
