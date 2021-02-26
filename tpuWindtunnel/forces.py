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

# ------------------------------------------------------------------------------
# Imported functions
# ------------------------------------------------------------------------------

import WindTunnelPostprocessing.aeroForces as aeroForces
import WindTunnelPostprocessing.modelProp as modelProp
import WindTunnelPostprocessing.scaling as scaling

import pyLEK.plotters.plot2D as plt

from pyLEK.helpers.pyExtras import getKeyList
from pyLEK.helpers.filemanager import delFilesInFolder

# ------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------     

def main():
    # Get name of input file
    # fname = sys.argv [1]
    fname = "C://Users//ac135564////bwSyncShare//MT Lais Windlastverteilung//21_Feb_Windkanaldaten//WindTunnelPostprocessing//T115_6_000.mat"

    # Clean up results folder
    # delFilesInFolder('T115_6/results')

    # Full scale building properties
    uH_f    = 37             # m/s       // Wind speed at z = H (50yr)
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

        # Time series of base forces
        t = np.linspace(0, buildProp.nT*buildProp.dT, buildProp.nT)
        F = buildAeroForces.BF_p
        plt.plot2D(t, F, "Time [s]", "Base Force " + dn + " [kN]", "Time series of base forces",
            ["dn"], showPlt=True)

if __name__ == '__main__':
    main()
    wait = input("Press Enter to continue.")

