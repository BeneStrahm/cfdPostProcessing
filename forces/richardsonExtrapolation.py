# -*- coding: utf-8 -*-
# ------------------------------------------------------------------------------
# Description:  Apply Richardson Extrapolation to forces to eval. grid convergence
# Author:       benedikt.strahm@ilek.uni-stuttgart.de
# Created:      2018-01-06
# Execution:    Import functions / collections (from cfdPostProcessing.forces import util)
#               or executing from command line (py filename.py)
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Sources
# ------------------------------------------------------------------------------
# Stern, Wilson, Coleman, Paterson - Comprehensive Approach to Verification and 
# Validation of CFD Simulations (2001), J. Fluids Eng. 

# ------------------------------------------------------------------------------
# Libraries
# ------------------------------------------------------------------------------
import numpy as np

import pyLEK.plotters.plot2D as plt
import os,sys,csv,h5py
from cfdPostProcessing.foamHelpers import readForces
from cfdPostProcessing.foamHelpers import readMoments

# ------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------


class convergenceVerification():
    # For a maximum of m = 3 solutions
    # Where m = 1 is the solution on the most fine parameters and 
    # m = 3 on the most coarse (e.g Grid, time stepping...)

    # See "Stern, Wilson, Coleman, Paterson - Comprehensive Approach
    # to Verification and Validation of CFD Simulations (2001)"
   
    def __init__(self): 
        

        # List of datapoints to be evaluated 
        self.convrgData     = [1,2,3]
        self.inputParam     = [1, 0.5, 0.25]
        self.delta_star_k_1 = "-"
        self.S_C            = "-"
        self.relError       = ["-", "-", "-"]

        # Calculate the refinement ratio Delta_x_k_n_+_1 / Delta_x_k_n
        self.r_k_21 = self.inputParam[1] / self.inputParam[0]
        self.r_k_32 = self.inputParam[2] / self.inputParam[1]

    def appendData(self, datapoint, m):
        # Insert datapoint at location m
        self.convrgData[m] = datapoint

    def evaluateConvergence(self):
        # Get the simulation result S_k, not corrected for iterative errors    
        S_k_1 = self.convrgData[0]
        S_k_2 = self.convrgData[1]
        S_k_3 = self.convrgData[2]

        # Calculate the solution changes epsilon_k
        self.epsilon_k_21 = S_k_2 - S_k_1
        self.epsilon_k_32 = S_k_3 - S_k_2

        # Solution change ratio
        self.R_k = (self.epsilon_k_21 / self.epsilon_k_32) / (self.r_k_21 / self.r_k_32)

        # Evaluate convergence condition
        if 0 < self.R_k < 1:
            self.convrgFlag = "Convergent"
        elif self.R_k <= 0:
            self.convrgFlag = "Oscillatory"
        elif self.R_k >= 1:
            self.convrgFlag = "Divergent"

    def estimateError(self):
        # Check convergence
        if self.convrgFlag == "Convergent":
            # Estimate order of accuracy rho_k
            # Initial guess 
            rho_k = np.log(self.epsilon_k_32 / self.epsilon_k_21 ) / \
                    np.log(np.mean([self.r_k_21, self.r_k_32]))

            # Iterative procedure 
            # rho_k_next = rho_k(i+1)
            # rho_k   = rho_k(i)
            # Stoping criteria: Relative Residual <= Value

            stopFlag = False

            while stopFlag == False:
                # Estimate next step
                rho_k_next = np.log(self.epsilon_k_32 / self.epsilon_k_21 ) /   \
                            np.log(self.r_k_21) + 1 / (np.log(self.r_k_21)) *     \
                            ( np.log(self.r_k_32 ** rho_k - 1) -                  \
                            np.log(self.r_k_21 ** rho_k - 1) )

                # Check criteria
                if abs((rho_k_next - rho_k) / rho_k_next) <= 0.001:
                    stopFlag = True

                # Update value with new iteration
                rho_k = rho_k_next

            # Assign to instance
            self.rho_k = rho_k
            
            # Estimated error
            self.delta_star_k_1 = self.epsilon_k_21 / (self.r_k_32 ** rho_k -1)

            # Corrected Solution results
            S_k_1 = self.convrgData[0]
            self.S_C = S_k_1 - self.delta_star_k_1

            # Relative Errors
            for key, S_k in enumerate(self.convrgData):
                S_C = self.S_C
                self.relError[key] = abs((S_C - S_k) / S_C)
    
    def writeHeader(self, outDir):
        # Open CSV-File in write mode
        outFile = outDir + "/ConvergenceTable.txt"
        
        with open(outFile, 'w') as txt:            
            txt_writer = csv.writer(txt, delimiter='\t', quotechar='"', quoting=csv.QUOTE_NONE)
            # Write refinement parameters
            txt_writer.writerow([
                "Grid A: Smallest grid cell size:",
                "Deltah3", 
                '{:4.3f}'.format(self.inputParam[2])
                ])
            txt_writer.writerow([
                "Grid B: Smallest grid cell size:",
                "Deltah2", 
                '{:4.3f}'.format(self.inputParam[1])
                ])
            txt_writer.writerow([
                "Grid C: Smallest grid cell size:",
                "Deltah1", 
                '{:4.3f}'.format(self.inputParam[0])
                ])
            txt_writer.writerow([
                ])
            txt_writer.writerow([
                "Refinement Ratio Grid A - B:",
                "rk32", 
                '{:4.3f}'.format(self.r_k_32)
                ])
            txt_writer.writerow([
                "Refinement Ratio Grid B - C:",
                "rk21", 
                '{:4.3f}'.format(self.r_k_21)
                ])
            txt_writer.writerow([
                ])

            # Write header of results
            # Characters per Tab: (['1___________14'])
            txt_writer.writerow([
                "Force component",                             
                "Parameter",  
                "Sk3: Result Grid A (Coarse)",   
                "Sk2: Result Grid B (Medium)",                                 
                "Sk1: Result Grid C (Fine)",
                "Rk: Convergence Ratio",
                "Convergence behaviour",
                "deltak1: Estimated Error of Sk1",
                "SC:Corrected Solution",
                "Relative Error Grid A",
                "Relative Error Grid B",
                "Relative Error Grid C"
                ])

    def writeConvergence(self, comp, designation, outDir, sT, eT):
        # Open CSV-File in append mode
        outFile = outDir + "/ConvergenceTable_" + str(int(sT)) + "_" + str(int(eT)) + ".txt"
        with open(outFile, 'a') as txt:            
            txt_writer = csv.writer(txt, delimiter='\t', quotechar='"', quoting=csv.QUOTE_NONE)

            # Check convergence
            # if True:
            if self.convrgFlag == "Convergent":
                # Write in results in row
                # Characters per Tab: (['1___________14'])

                txt_writer.writerow([str(comp),                             
                    str(designation),                                       
                    '{: 10.3f}'.format(self.convrgData[2]),                 
                    '{: 10.3f}'.format(self.convrgData[1]),                 
                    '{: 10.3f}'.format(self.convrgData[0]),                 
                    '{: 6.3f}'.format(self.R_k),
                    str(self.convrgFlag),
                    '{: 8.3f}'.format(self.delta_star_k_1),
                    '{: 10.3f}'.format(self.S_C),
                    '{: 6.1f}'.format(self.relError[2]*100),
                    '{: 6.1f}'.format(self.relError[1]*100),
                    '{: 6.1f}'.format(self.relError[0]*100),
                    str(designation),                                       
                    '{: 10.3f}'.format(self.convrgData[2]),                 
                    '{: 10.3f}'.format(self.convrgData[1]),                 
                    '{: 10.3f}'.format(self.convrgData[0]),                 
                    '{: 6.3f}'.format(self.R_k),
                    str(self.convrgFlag),
                    self.delta_star_k_1,
                    self.S_C,
                    self.relError[2],
                    self.relError[1],
                    self.relError[0]
                    ])

    def plotConvergence(self, outDir, sT, eT):

        dir_fileName = outDir + "Baseshear"
        x = self.inputParam
        y = self.convrgData

        if self.convrgFlag == "Convergent":
            hLine = self.S_C
            hText = 'estimate'
        else:
            hLine = None
            hText = None

        # # 1) For Mean
        # if "Mean" in designation:
        #     if "F" in comp:
        #         ylabel  = r"\$\overbar{F_{" + r"{}".format(comp[1]) + r"}}~[MN]\$" 
        #         title   = r"Convergence Mean F \textsubscript{" + r"{}".format(comp)[1] + r"}"
        #     elif "M" in comp:
        #         ylabel  = r"\$\overbar{M_{" + r"{}".format(comp[1]) + r"}}~[MNm]\$" 
        #         title   = r"Convergence Mean M \textsubscript{" + r"{}".format(comp)[1] + r"}"
        
        # # 2) For RMS
        # elif "RMS" in designation:
        #     if "F" in comp:
        #         ylabel  = r"\$F_{" + r"{}".format(comp[1]) + r",rms}~[MN]\$" 
        #         title   = r"Convergence RMS F \textsubscript{" + r"{}".format(comp)[1] + r"}"
        #     elif "M" in comp:
        #         ylabel  = r"\$M_{" + r"{}".format(comp[1]) + r",rms}~[MNm]\$" 
        #         title   = r"Convergence RMS M \textsubscript{" + r"{}".format(comp)[1] + r"}"
        
        # # 3) For Std
        # elif "Std" in designation:
        #     if "F" in comp:
        #         ylabel  = r"\$F_{" + r"{}".format(comp[1]) + r",}~[MN]\$" 
        #         title   = r"Convergence Std F \textsubscript{" + r"{}".format(comp)[1] + r"}"
            # elif "M" in comp:
            #     ylabel  = r"\$M_{" + r"{}".format(comp[1]) + r",std}~[MNm]\$" 
            #     title   = r"Convergence Std M \textsubscript{" + r"{}".format(comp)[1] + r"}"
        
        # # 4) For Largest20Values
        # elif "Largest20" in designation:
        #     if "F" in comp:
        #         ylabel  = r"\$F_{" + r"{}".format(comp[1]) + r",max}~[MN]\$" 
        #         title   = r"Convergence Max F \textsubscript{" + r"{}".format(comp)[1] + r"}"
        #     elif "M" in comp:
        #         ylabel  = r"\$M_{" + r"{}".format(comp[1]) + r",max}~[MNm]\$" 
        #         title   = r"Convergence Max M \textsubscript{" + r"{}".format(comp)[1] + r"}"
        
        # # 5) For Smallest20Values
        # elif "Smallest20" in designation:
        #     if "F" in comp:
        #         ylabel  = r"\$F_{" + r"{}".format(comp[1]) + r",min}~[MN]\$" 
        #         title   = r"Convergence Min F \textsubscript{" + r"{}".format(comp)[1] + r"}"
        #     elif "M" in comp:
        #         ylabel  = r"\$M_{" + r"{}".format(comp[1]) + r",min}~[MNm]\$" 
        #         title   = r"Convergence Min M \textsubscript{" + r"{}".format(comp)[1] + r"}"



        

        xlim = [max(x)+0.1*max(x),min(x)-0.5*min(x)]
        ylim = [0, 0.5]
        # Change to current file location
        os.chdir(os.path.dirname(sys.argv[0]))

        style_dict = {"lines.markersize": 8,"savefig.format": "svg"}
        xlabel = 'Zellgröße'
        ylabel  = r"Fx Std [MN]" 
        title   = r"Convergence Std"

        plt.plot2D(x, y, xlabel=xlabel, ylabel=ylabel, title=title, 
               dir_fileName=dir_fileName, xlim=xlim, ylim=ylim,
               style_dict=style_dict,
               variation='marker',
               savePlt=True, showPlt=True)

    
# ------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------

def main():
    # Get ouput directory
    outDir = '/media/dani/linuxHDD/openfoam/simpleFoam/testing/'

    # Collect the components
    fname0 = '/media/dani/linuxHDD/openfoam/simpleFoam/testing/1_conv_ref0/postProcessing/forces/0/force.dat'
    fname1 = '/media/dani/linuxHDD/openfoam/simpleFoam/testing/1_conv_ref1/postProcessing/forces/0/force.dat'
    fname2 = '/media/dani/linuxHDD/openfoam/simpleFoam/testing/1_conv_ref2/postProcessing/forces/0/force.dat'

    interpForces0 = readForces.importForces(fname0)
    interpForces1 = readForces.importForces(fname1)
    interpForces2 = readForces.importForces(fname2)
    interp = [interpForces0, interpForces1, interpForces2]

    # Specify Cut-Off time
    print("Specify time range to plot")
    sCutT = 200
    eCutT = 250

    # Set-Up CSV-Table
    tableHeader = convergenceVerification()
    tableHeader.writeHeader(outDir)
    
    # Loop through components
    # meanBF  = convergenceVerification()
    # rmsBF   = convergenceVerification()
    stdBF   = convergenceVerification()
    # maxBF   = convergenceVerification()
    # minBF   = convergenceVerification()

    m = 0
    while m <= 2:
            # Get data to plot
            dT = interp[m][0][1]-interp[m][0][0]
            BF = interp[m][2]/ (10 ** 6)    
            
            # Cut the array to desired range 
            l = BF[int(sCutT/dT):int(eCutT/dT)]
            # meanBF.appendData(np.mean(l), m)
            # rmsBF.appendData(np.sqrt(np.mean(np.square(l))), m)
            stdBF.appendData(np.std(l), m)  
            print (np.std(l))
            # maxBF.appendData(np.average(hq.nlargest(20, l)), m) 
            # minBF.appendData(np.average(hq.nsmallest(20, l)), m) 
            m +=1


    # Evaluate convergence(s)

    # meanBF.evaluateConvergence()
    # rmsBF.evaluateConvergence()
    stdBF.evaluateConvergence()
    # maxBF.evaluateConvergence()
    # minBF.evaluateConvergence()

    # Estimate the error(s)

    # meanBF.estimateError()
    # rmsBF.estimateError()
    stdBF.estimateError()
    # maxBF.estimateError()
    # minBF.estimateError()

    # Write convergence to txt-file

    # meanBF.writeConvergence(comp,"Mean", outDir, sCutT, eCutT)
    # rmsBF.writeConvergence(comp,"RMS", outDir, sCutT, eCutT)
    # stdBF.writeConvergence(comp,"Std", outDir, sCutT, eCutT)
    # maxBF.writeConvergence(comp,"Largest20", outDir, sCutT, eCutT)
    # minBF.writeConvergence(comp,"Smallest20", outDir, sCutT, eCutT)

    # Plot convergence (+ estimated error)
    # meanBF.plotConvergence(comp, "Mean", outDir, sCutT, eCutT)
    # rmsBF.plotConvergence(comp, "RMS", outDir, sCutT, eCutT)
    stdBF.plotConvergence(outDir, sCutT, eCutT)
    # maxBF.plotConvergence(comp,"Largest20", outDir, sCutT, eCutT)
    # minBF.plotConvergence(comp,"Smallest20", outDir, sCutT, eCutT)

if __name__ == '__main__':
    main()
