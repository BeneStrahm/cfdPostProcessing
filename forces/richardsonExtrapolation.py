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
# Celik - Procedure for Estimation and Reporting of  Discretization Error in 
# CFD Applications (2008), J. Fluids Eng. 

# ------------------------------------------------------------------------------
# Libraries
# ------------------------------------------------------------------------------
from matplotlib.pyplot import hlines
import numpy as np

import pyLEK.plotters.plot2D as plt
import os,sys,csv,h5py
from cfdPostProcessing.foamHelpers import readForces
from cfdPostProcessing.foamHelpers import readMoments
import scipy.optimize as optimize

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
        self.S_k_m          = [1,2,3]
        self.delta_x_k_m    = [1,2,3]
        self.delta_star_k_1 = "-"
        self.S_C            = "-"
        self.relError       = ["-", "-", "-"]

    def appendData(self, S_k_m, delta_x_k_m, m):
        # Insert datapoint at location m
        self.S_k_m[m] = S_k_m
        self.delta_x_k_m[m] = delta_x_k_m

    def evaluateConvergence(self):
        # Get the simulation result S_k, not corrected for iterative errors 
        self.delta_x_k_1 = self.delta_x_k_m[0]
        self.delta_x_k_2 = self.delta_x_k_m[1]
        self.delta_x_k_3 = self.delta_x_k_m[2]

        # Calculate the refinement ratio Delta_x_k_n_+_1 / Delta_x_k_n
        self.r_k_21 = self.delta_x_k_2 / self.delta_x_k_1
        self.r_k_32 = self.delta_x_k_3 / self.delta_x_k_2

        # Get the simulation result S_k, not corrected for iterative errors    
        self.S_k_1 = self.S_k_m[0]
        self.S_k_2 = self.S_k_m[1]
        self.S_k_3 = self.S_k_m[2]

        # Calculate the solution changes epsilon_k
        self.epsilon_k_21 = self.S_k_2 - self.S_k_1
        self.epsilon_k_32 = self.S_k_3 - self.S_k_2

        # Relative errors
        self.e_a_21 = abs((self.epsilon_k_21) / self.S_k_1)
        self.e_a_32 = abs((self.epsilon_k_32) / self.S_k_2)

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
                    np.log(self.r_k_21)
            s     = np.sign(self.epsilon_k_32/self.epsilon_k_21)

            # Iterative procedure 
            # rho_k_next = rho_k(i+1)
            # rho_k   = rho_k(i)
            # Stoping criteria: Relative Residual <= Value
            stopFlag = False

            while stopFlag == False:
                # Estimate next step (Improved meth., see Celik)
                rho_k_next = 1 / np.log(self.r_k_21) * \
                             np.abs(np.log(
                                    np.abs(self.epsilon_k_32/self.epsilon_k_21)
                                )) \
                            + np.log(
                                (self.r_k_21**rho_k - s) / (self.r_k_32**rho_k - s)  
                            )
                            

                # Check criteria
                if abs((rho_k_next - rho_k) / rho_k_next) <= 0.001:
                    stopFlag = True

                # Update value with new iteration
                rho_k = rho_k_next

            # Assign to instance
            self.rho_k = rho_k
            
            # Estimated error 
            self.delta_star_k_1 = self.epsilon_k_21 / (self.r_k_32 ** rho_k -1)

            # Grid Convergence Index (GCI)
            self.F_S = 1.25      # Safety Factor (Fs=1.25 for comparisons over three or more grids)
            self.GCI_21 = self.F_S * self.e_a_21 / (self.r_k_21 ** rho_k -1)

            # Corrected Solution results
            self.S_k_1 = self.S_k_m[0]
            self.S_C = self.S_k_1 - self.delta_star_k_1

            # Extrapolated relative errors
            for key, S_k in enumerate(self.S_k_m):
                S_C = self.S_C
                self.relError[key] = abs((S_C - S_k) / S_C)

    def writeConvergence(self, comp, designation, outDir):
        # Open CSV-File in append mode
        outFile = outDir + "/ConvergenceTable.txt"
        with open(outFile, 'a') as txt:            
            txt_writer = csv.writer(txt, delimiter='\t', quotechar='"', quoting=csv.QUOTE_NONE)
            # Write refinement parameters
            txt_writer.writerow([
                "Grid A: Smallest grid cell size:",
                "Deltah3", 
                '{:4.3f}'.format(self.delta_x_k_m[2])
                ])
            txt_writer.writerow([
                "Grid B: Smallest grid cell size:",
                "Deltah2", 
                '{:4.3f}'.format(self.delta_x_k_m[1])
                ])
            txt_writer.writerow([
                "Grid C: Smallest grid cell size:",
                "Deltah1", 
                '{:4.3f}'.format(self.delta_x_k_m[0])
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

            # Check convergence
            # if True:
            if self.convrgFlag == "Convergent":
                # Write in results in row
                # Characters per Tab: (['1___________14'])

                txt_writer.writerow([str(comp),                             
                    str(designation),                                       
                    '{: 10.3f}'.format(self.S_k_m[2]),                 
                    '{: 10.3f}'.format(self.S_k_m[1]),                 
                    '{: 10.3f}'.format(self.S_k_m[0]),                 
                    '{: 6.3f}'.format(self.R_k),
                    str(self.convrgFlag),
                    '{: 8.3f}'.format(self.delta_star_k_1),
                    '{: 10.3f}'.format(self.S_C),
                    '{: 6.1f}'.format(self.relError[2]*100),
                    '{: 6.1f}'.format(self.relError[1]*100),
                    '{: 6.1f}'.format(self.relError[0]*100),
                    str(designation),                                       
                    '{: 10.3f}'.format(self.S_k_m[2]),                 
                    '{: 10.3f}'.format(self.S_k_m[1]),                 
                    '{: 10.3f}'.format(self.S_k_m[0]),                 
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
        x = self.delta_x_k_m
        y = self.S_k_m

        if self.convrgFlag == "Convergent":
            hLines = [self.S_C]
            hTexts = ['estimate']
        else:
            hLines = None
            hTexts = None

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



        

        xlim = [max(x)+max(x)*0.2,0]
        ylim = [-0.2,max(y)+max(y)*0.5]
        # Change to current file location
        os.chdir(os.path.dirname(sys.argv[0]))

        style_dict = {"lines.linewidth": 0,"lines.markersize": 8,"savefig.format": "svg","lines.marker":"x"}
        xlabel = 'Zellanzahl'
        ylabel  = r"Fx Std [MN]" 
        title   = r"Convergence Std"
        

        plt.plot2D(x, y, xlabel=xlabel, ylabel=ylabel, title=title, 
               dir_fileName=dir_fileName, xlim=xlim,# ylim=ylim,
               hLines=hLines, hTexts=hTexts,
               style_dict=style_dict,
               variation='color',
               savePlt=True, showPlt=True)

    
# ------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------

def main():
    # Get ouput directory
    outDir = '/media/dani/linuxHDD/openfoam/simpleFoam/testing/'
    
    #bs
    # outDir = 'C:/Users/bstra/GitHub/cfdPostProcessing/forces/sampleData/'

    # Step 1: Define a representative grid size h
    # -------
    # Where m = 1  most fine parameters and m = 3 on the most coarse 
    no_of_cells     = [4924000, 3207000, 1000000]
    refinement_vol  = 3840 * 1600 * 960
    grid_size_h     = [(refinement_vol / x) ** (1/3) for x in no_of_cells]

    # Collect the components
    fname0 = '/media/dani/linuxHDD/openfoam/simpleFoam/testing/1_conv_ref0/postProcessing/forces/0/force.dat'
    fname1 = '/media/dani/linuxHDD/openfoam/simpleFoam/testing/1_conv_ref1/postProcessing/forces/0/force.dat'
    fname2 = '/media/dani/linuxHDD/openfoam/simpleFoam/testing/1_conv_ref2/postProcessing/forces/0/force.dat'

    # bs:
    # fname0 = 'C:/Users/bstra/GitHub/cfdPostProcessing/forces/sampleData/force_ref0.dat'
    # fname1 = 'C:/Users/bstra/GitHub/cfdPostProcessing/forces/sampleData/force_ref1.dat'
    # fname2 = 'C:/Users/bstra/GitHub/cfdPostProcessing/forces/sampleData/force_ref2.dat'

    interpForces2 = readForces.importForces(fname0)
    interpForces1 = readForces.importForces(fname1)
    interpForces0 = readForces.importForces(fname2)
    interp = [interpForces0, interpForces1, interpForces2]

    # Specify Cut-Off time
    # print("Specify time range to plot")
    sCutT = 200
    eCutT = 250
 
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
            BF = interp[m][1]/ (10 ** 6)    
            
            # Cut the array to desired range 
            l = BF[int(sCutT/dT):int(eCutT/dT)]
            # meanBF.appendData(np.mean(l), m)
            # rmsBF.appendData(np.sqrt(np.mean(np.square(l))), m)
            stdBF.appendData(np.std(l), grid_size_h[m], m)   
            # print (np.sqrt(np.mean(np.square(l))))
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
    stdBF.writeConvergence("Fx","Std", outDir)
    # maxBF.writeConvergence(comp,"Largest20", outDir, sCutT, eCutT)
    # minBF.writeConvergence(comp,"Smallest20", outDir, sCutT, eCutT)

    # Plot convergence (+ estimated error)
    # meanBF.plotConvergence(outDir, sCutT, eCutT)
    # rmsBF.plotConvergence(outDir, sCutT, eCutT)
    stdBF.plotConvergence(outDir, sCutT, eCutT)
    # maxBF.plotConvergence(comp,"Largest20", outDir, sCutT, eCutT)
    # minBF.plotConvergence(comp,"Smallest20", outDir, sCutT, eCutT)

if __name__ == '__main__':
    main()
