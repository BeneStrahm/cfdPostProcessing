# ------------------------------------------------------------------------------
# Project       : Master Thesis
# Description   : Evaluate & Plot simulation run time (for AMR Simulations)
# Author        : benedikt.strahm@tum.de
# Created       : 2021-04-21
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Libraries
# ------------------------------------------------------------------------------

import os
import sys
import re
import csv
import numpy as np

# ------------------------------------------------------------------------------
# Imported functions
# ------------------------------------------------------------------------------

import pyLEK.plotters.plot2D as plt

# ------------------------------------------------------------------------------
# Classes
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------


def safeReturnData(data):
    try:
        i = float(data)
    except ValueError:
        i = None
    return i


def main():
    # Get log file
    import tkinter as tk
    from tkinter import filedialog

    # root = tk.Tk()
    # root.withdraw()

    # inDir = filedialog.askdirectory()
    inDir = 'C:/Users/ac135564/GitHub/cfdPostProcessing/runtime/sampleData'

    # outDir      = createOurDir.main("../Results/OpenFOAM/" + "RunTime/" + dObject.Name)
    # outFile     = outDir + "/runTimeEval.txt"
    inDirRun = inDir + "/logRun"
    inFileMesh = inDir + "/logMesh/checkMesh.log"

    # Get list of files to be imported
    logFiles = []
    for logFile in os.listdir(inDirRun):
        if "Foam_Run" in logFile:
            logFiles.append(logFile)

    # Sort files
    logFiles.sort()

    # Initialize (T)otal run time
    T_Sim = []
    T_CPU = []
    T_Clock = []
    deltaT = []

    # Initialize CFL numbers
    CFL_max = []
    CFL_mean = []

    # Initialize refinement levels
    CFL_max = []
    CFL_mean = []

    # Initialize number of cells
    n_cells_ref = []
    T_cells_ref = []
    n_cells_unref = []
    T_cells_unref = []

    # Reading from files
    # ------------------

    # Loop over files
    for logFile in logFiles:

        # Set initial restart flag when new logFile is opened
        flagRestart = True

        path = os.path.join(inDirRun, logFile)
        with open(os.path.join(path), 'r') as dat:
            # Loop over lines
            for line in dat:

                # Check if keyword is in line
                if "Time = " in line:

                    # Split lines at " " and replace newline signs
                    line = line.replace("\n", "")
                    line = line.split(" ")

                    # Get simulation-, cpu- and clocktimes
                    if line[0] == "Time":

                        # Write current time step
                        T_Sim.append(safeReturnData(line[2]))

                    elif line[0] == "ExecutionTime":

                        # Write CPU / Clock times
                        T_CPU.append(safeReturnData(line[2]))
                        T_Clock.append(safeReturnData(line[7]))

                elif "deltaT = " in line:

                    # Split lines at " " and replace newline signs
                    line = line.replace("\n", "")
                    line = line.split(" ")

                    # Get current deltaT
                    if line[0] == "deltaT":

                        # Write current deltaT
                        deltaT.append(safeReturnData(line[2]))

                elif "Courant Number" in line:
                    # Split lines at " " and replace newline signs
                    line = line.replace("\n", "")
                    line = line.split(" ")

                    # Write Courant (CFL) Number
                    CFL_max.append(safeReturnData(line[5]))
                    CFL_mean.append(safeReturnData(line[3]))

                elif "Refined from" in line:
                    # Split lines at " " and replace newline signs
                    line = line.replace("\n", "")
                    line = line.split(" ")

                    # Courant (CFL) Number to
                    n_cells_ref.append(safeReturnData(line[4]))
                    T_cells_ref.append(T_Sim[-1])

                elif "Unrefined from" in line:
                    # Split lines at " " and replace newline signs
                    line = line.replace("\n", "")
                    line = line.split(" ")

                    # Courant (CFL) Number to
                    n_cells_unref.append(safeReturnData(line[4]))
                    T_cells_unref.append(T_Sim[-1])

    # Post processing
    # ------------------
    # Calculate run time per time step dT_CPU/dT
    # Note that this approximated "derivative" has size n-1 where n is your array/list size.
    dT_CPU_dT = np.diff(T_CPU)

    # # Write results
    # with open(outFile, 'w') as txt:
    #     txt_writer = csv.writer(txt, delimiter='\t', quotechar='"', quoting=csv.QUOTE_NONE)
    #     txt_writer.writerow(["Simulation Time:",                np.max(T_Sim)])
    #     txt_writer.writerow(["# of time steps:",                np.max(T_Sim)])
    #     txt_writer.writerow(["Delta Time:",                     np.max(T_Sim)])
    #     txt_writer.writerow(["Real CPU Time:",                  np.max(T_CPU)])
    #     txt_writer.writerow(["Real Clock Time:",                np.max(T_Clock)])
    #     txt_writer.writerow(["Avg CPU Time / time step:",       np.max(T_CPU)/nT])
    #     txt_writer.writerow(["Avg CPU Time / time step / cell:",np.max(T_CPU)/nT/numOfCells])
    #     txt_writer.writerow([])
    #     txt_writer.writerow(["Run Time between 200 - 500 steps"])
    #     txt_writer.writerow(["CPU Time:"                       ,T_CPU[499] - T_CPU[199]])
    #     txt_writer.writerow(["Time Steps:"                     ,(T_Sim[499] - T_Sim[199])/dT])
    #     txt_writer.writerow(["Avg CPU Time / time step:"       ,(T_CPU[499] - T_CPU[199])/(T_Sim[499] - T_Sim[199]) * dT])

    # Plot results
    # ------------------
    # x = T_Sim
    # y = [T_CPU, T_Clock]

    # xlabel = "\Wind Tunnel Time [s]\\"
    # ylabel = "\Computational Time [s]\\"
    # title = "Computational Run Time Evaluation"
    # legend = ["CPU Time", "Clock Time"]

    # # Change to current file location
    # os.chdir(os.path.dirname(sys.argv[0]))

    # dir_fileName = "sampledata/computational_run_time"

    # figSize = plt.plotHelpers.calcFigSize()

    # style_dict = {"lines.linewidth": 0.7, "figure.figsize": figSize}

    # plt.plot2D(x, y, xlabel=xlabel, ylabel=ylabel, title=title, legend=legend,
    #            style_dict=style_dict, dir_fileName=dir_fileName,
    #            savePlt=False, showPlt=True)

    x = T_Sim[-1]
    y = [dT_CPU_dT]

    xlabel = "\Wind Tunnel Time [s]\\"
    ylabel = "\Computational Time [s]\\"
    title = "Computational Run Time Evaluation"
    legend = ["CPU Time"]

    # Change to current file location
    os.chdir(os.path.dirname(sys.argv[0]))

    dir_fileName = "sampledata/computational_run_time"

    figSize = plt.plotHelpers.calcFigSize()

    style_dict = {"lines.linewidth": 0.7, "figure.figsize": figSize}

    plt.plot2D(x, y, xlabel=xlabel, ylabel=ylabel, title=title, legend=legend,
               style_dict=style_dict, dir_fileName=dir_fileName, variation='marker',
               savePlt=False, showPlt=True)

    # x = [T_Sim]
    # y = [CFL_mean, CFL_max]

    # xlabel = "\Wind Tunnel Time [s]\\"
    # ylabel = "\CFL\\"
    # title = "Evolution of Courant-Friedrich-Lewis Number"
    # dir_fileName = outDir + "CFL"
    # legend = ["mean CFL", "max CFL"]

    # xlim = []
    # ylim = []

    # style_dict = {"figure.figsize":"5.83, 3.28" }

    # plt.plot2D(x, y, xlabel, ylabel, title, legend, dir_fileName,
    #             xlim=xlim, ylim=ylim, xscale='linear', yscale='linear',
    #             style_dict=style_dict, colorScheme='TUM', variation='color',
    #             transparent=False, alphaLines=0.7)


if __name__ == '__main__':
    main()
