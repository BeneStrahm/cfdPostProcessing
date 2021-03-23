#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 17:53:18 2021

@author: iagaxtma
"""

#%%
###############################################################################
#Modules
###############################################################################

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
import sys
import re
import os
import glob
import pyvista as pv
import pandas as pd


#%%
###############################################################################
#Plotting stuff SECTION
###############################################################################

def latexify():
  scale=1
  page_width_pt=485*scale #pt
  inches_per_pt=1/72.27
  golden_mean = (np.sqrt(5.0)-1.0)/2.0           #was /2.0 Aesthetic ratio (you could change this)
  fig_width=page_width_pt*inches_per_pt
  fig_height=fig_width*golden_mean  
  fig_size=[fig_width,fig_height]

  params = {'backend': 'ps',
            #'text.latex.preamble': [r'\usepackage{serif}'],
            'axes.labelsize': 10*scale, # fontsize for x and y labels (was 10)
            'axes.titlesize': 10*scale,
            #'text.fontsize': 10*scale, # was 10
            'legend.fontsize': 10*scale, # was 10
            'xtick.labelsize': 10*scale,
            'ytick.labelsize': 10*scale,
            'text.usetex': True,
            'figure.figsize': fig_size,
            'axes.linewidth':0.5,
            'xtick.direction': 'out',
            'ytick.direction': 'out',
            'font.family': 'serif',
            'grid.color'        : 'lightgray',      # grid color, 'lightgrey','lightgray','silver'
            'grid.linestyle'    : '-',
            'grid.linewidth'    : 0.5,      # in points
            'grid.alpha'        : 1.0,  
            "xtick.minor.visible"  : False# transparency, between 0.0 and 1.0
  }
  plt.rcParams.update(params)
   
def savefig(filename):
    #plt.savefig('{}.pgf'.format(filename),dpi=1000,bbox_inches='tight')
    plt.savefig('{}.pdf'.format(filename),dpi=600,bbox_inches='tight')
    #tikz_save('{}.tikz'.format(filename),
           #figureheight = '\\figureheight',
           #figurewidth = '\\figurewidth')
           
def savefigPNG(filename):
    plt.savefig('{}.png'.format(filename),dpi=300,bbox_inches='tight')
    

new_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728','#9467bd', '#8c564b', '#e377c2', '#7f7f7f','#bcbd22', '#17becf']

#%%
###############################################################################
#Defintiion Functions
###############################################################################

def list_files(dir):
    r = []
    for root, dirs, files in os.walk(dir):
        for name in files:
            if name.endswith("vtp"):
                r.append(os.path.join(root, name))
    r.sort(key=os.path.getctime)
    return r


def loadData(filename,fieldScalar,ofDisc=False)    :
    print('\t \t Loading and Interpolating file:' + str(filename) + '\t field: ' + str(fieldScalar))
    data    = pv.read(filename)
    mesh=data.points
    x,y,z=mesh[:,0],mesh[:,1],mesh[:,2]

    p=data.point_arrays.get_array(fieldScalar)

    nx,ny,nz=d*2+1,d*2+1,h+1
       
    if (ofDisc==False):     
        xs1D = np.linspace(x.min(),x.max(),nx,endpoint=True)
        ys1D = np.linspace(y.min(),y.max(),ny,endpoint=True)
        zs1D = np.linspace(z.min(),z.max(),nz,endpoint=True)
        print('\t \t using user discretization nx,ny,nz = ' + str(nx) + ',' + str(ny) + ',' + str(nz))
    
    x0=np.array([x.min()])[0]
    x1=np.array([x.max()])[0]
    y0=np.array([y.min()])[0]
    y1=np.array([y.max()])[0]
    z0=np.array([z.min()])[0]
    z1=np.array([z.max()])[0]
    
    
    print('\t \t interpolating and creating linear cube map')
    [xs,ys,zs] = np.meshgrid(np.array(nx*[x0]),ys1D,zs1D)
    piFront=interpolate.griddata((x,y,z),p,(xs,ys,zs),method='nearest')[:,0,:]
    [xs,ys,zs] = np.meshgrid(np.array(nx*[x1]),ys1D,zs1D)
    piBack=-1*interpolate.griddata((x,y,z),p,(xs,ys,zs),method='nearest')[:,0,:]
    [xs,ys,zs] = np.meshgrid(xs1D,np.array(ny*[y0]),zs1D)
    piSide1=interpolate.griddata((x,y,z),p,(xs,ys,zs),method='nearest')[0,:,:]
    # piSide1=np.flip(piSide1,axis=0)
    [xs,ys,zs] = np.meshgrid(xs1D,np.array(ny*[y1]),zs1D)
    piSide2=-1*interpolate.griddata((x,y,z),p,(xs,ys,zs),method='nearest')[0,:,:]
    
    pFB=np.concatenate((piFront,piBack))
    pSS=np.concatenate((piSide1,piSide2))
    pAll=np.concatenate((piSide1,piFront,piSide2,piBack))
    
    xFB=np.linspace(0.5*y0,y1+0.5*y1,pFB.shape[0])
    xSS=np.linspace(y0,y1,pSS.shape[0])
    xAll=np.linspace(y0-0.5*y1,2*y1+0.5*y1,pAll.shape[0])
    zAll=np.linspace(z0,z1,pAll.shape[1])
  
    pAll=np.transpose(pAll)
    pFB=np.transpose(pFB)
    pSS=np.transpose(pSS)
     
    print('\t ... done')
    return xFB,xSS,xAll,zAll,pFB,pSS,pAll

def pressureMaxOverHeight(x,y):
    zmaxH=x
    pmaxH=np.max(y,axis=1)
    return zmaxH,pmaxH


def FFT(write):
    
    ###################################################
    # Temporal FFT for the real part and the amplitude of the spatial FFT
    # For the wave number from 0 to Kfft
    ###################################################

    print ("starting FFT for wave No from 0 to "+str(Kfft))
    """ By far the most efficient method """
    
    start = time()
    
    import multiprocessing
    nthread = multiprocessing.cpu_count()

    A=np.zeros((fftAmp.shape[0],fftAmp.shape[1],fftAmp.shape[2],Kfft+1),dtype='complex')    
    AmpTfftAmp=np.zeros((int(fftAmp.shape[0]/2),fftAmp.shape[1],fftAmp.shape[2],Kfft+1),dtype='complex')
    AmpTfftReal=np.zeros((int(fftAmp.shape[0]/2),fftAmp.shape[1],fftAmp.shape[2],Kfft+1),dtype='complex')
            
    RealTfftAmp=np.zeros((int(fftAmp.shape[0]/2),fftAmp.shape[1],fftAmp.shape[2],Kfft+1),dtype='complex')
    RealTfftReal=np.zeros((int(fftAmp.shape[0]/2),fftAmp.shape[1],fftAmp.shape[2],Kfft+1),dtype='complex')

    import pyfftw

    A = pyfftw.interfaces.numpy_fft.fft(fftAmp[:,:,:,0:Kfft+1], axis=0, threads=nthread)
    AmpTfftAmp = np.abs(A[0:int(fftAmp.shape[0]/2),:,:,:])
    AmpTfftReal = np.real(A[0:int(fftAmp.shape[0]/2),:,:,:])
    
    A = pyfftw.interfaces.numpy_fft.fft(fftReal[:,:,:,0:Kfft+1], axis=0, threads=nthread)
    RealTfftAmp = np.abs(A[0:int(fftReal.shape[0]/2),:,:,:])
    RealTfftReal = np.real(A[0:int(fftReal.shape[0]/2),:,:,:])
   
    end = time()
    print("Time for FFT ("+str(Kfft+1)+" wave numbers) "+ str(end - start))
    
    ###################################################
     #Saving the results of the FFT to OD
    ###################################################    

    if write:
        np.save(OD+'AmpTfftAmp',AmpTfftAmp)
        np.save(OD+'AmpTfftReal',AmpTfftReal)
        np.save(OD+'RealTfftAmp',RealTfftAmp)
        np.save(OD+'RealTfftReal',RealTfftReal)
        
    print('data format: Tfft[freq,x,y,waveNo]')
    return AmpTfftAmp, AmpTfftReal, RealTfftAmp, RealTfftReal

    
def chunks(l, n):
    n = max(1, n)
    return (l[i:i+n] for i in range(0, len(l), n))

def createFloorData(heightFloor,z,data):
    nFloors=round(z.max()/heightFloor) #in [m]
    chunkFloors=np.array_split(z, nFloors)
    chunkFloorIndex=[]
    for i in range(0,len(chunkFloors)):
        index=[]
        for j in range(0,len(chunkFloors[i])):
            index.append(np.argwhere(z==chunkFloors[i][j]))
        chunkFloorIndex.append(index)
    
    chunkFloors_z0=[chunkFloors[i][0] for i in range(0,len(chunkFloors))]
    chunkFloors_z1=[chunkFloors[i][-1] for i in range(0,len(chunkFloors))]
    chunkFloors_index0=[chunkFloorIndex[i][0] for i in range(0,len(chunkFloorIndex))]
    chunkFloors_index1=[chunkFloorIndex[i][-1] for i in range(0,len(chunkFloorIndex))]
    
    dataFloors=[]
    dataFloorsMean=[]
    for i in range(0,len(pFB)):
        temp=[]
        temp2=[]
        for j in range(0,int(nFloors)):
            temp.append(data[i][int(chunkFloors_index0[j]):int(chunkFloors_index1[j]),:])
            temp2.append(np.mean(data[i][int(chunkFloors_index0[j]):int(chunkFloors_index1[j]),:]))
        dataFloors.append(temp)
        dataFloorsMean.append(temp2)
    
    dataFloorsMean=np.transpose(np.array(dataFloorsMean))
    dataFloorsMean=dataFloorsMean*np.diff(z)[0]*np.diff(xFB)[0]*rho_inf
    zFloors=np.arange(0,nFloors,1)
    return zFloors, dataFloorsMean



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




#%%

###############################################################################
#INPUT SECTION
###############################################################################

#Control
interpolateData=True
loadinterpolatedData=True
createInstantAnimation=False
#Data Input
caseName='1_conv_ref1' 
rootPath='/media/dani/linuxHDD/openfoam/simpleFoam/testing/postProcessing/raw/' + caseName
filename='building_wall.vtp'                        #name of data file after inpolating data=160                                               #height
d=32                                                #length

#Global time data
t0,t1,dt=60,250,1
time=np.arange(t0,t1+dt,dt)
timeZeroed=time-t0                                  #Zero time

#Ambient stuff
rho_inf=1.18                                        #rho_inf
p_inf=101325                                        #pressure inf

heightFloor=4                                        #height floor [m]


#%%
###############################################################################
#Preprocessing data
###############################################################################

fileList=list_files(rootPath)
fields=[]

if (interpolateData==True):
#Create FileList Dictionary
    for i in range(0,len(fileList)):
        print('Preprocessing file: ' + str(i) + '/' + str(len(fileList)) + '\t time: ' + str(time[i]))
        xFB,xSS,xAll,zAll,pFB,pSS,pAll=loadData(fileList[i],'p',ofDisc=False)
        fields.append([timeZeroed[i],xFB,xSS,xAll,zAll,pFB,pSS,pAll])
        del xFB,xSS,xAll,zAll,pFB,pSS,pAll
        
    print('Converting and saving all field data as pickle file:')        
    import pandas as pd
    df = pd.DataFrame(fields)
    df.to_pickle(caseName)
    
if (loadinterpolatedData==True): 
    df=pd.read_pickle(caseName)

xFB=    df[1][0]
xSS=    df[2][0]
xAll=   df[3][0]
z=      df[4][0]
pFB=    df[5]*rho_inf
pSS=    df[6]*rho_inf
pALL=   df[7]*rho_inf
p = pALL
x = xAll
#%%Plotting contour data and saving
animationPath='./instantField/'
if not os.path.exists(animationPath):
    os.makedirs(animationPath)

if (createInstantAnimation==True):   
    levels_p=np.linspace(-1000,1000,64)
    for i in range(0,len(p)):
        print('Creating Animation of instant field ' + str(i) + '/' + str(len(fileList)))
        latexify()
        fig1, (ax1) =  plt.subplots(1)
        axlist=[ax1]
        surf=ax1.contourf(x,z,p[i]-p_inf,levels=levels_p,extend='both',cmap='coolwarm')
        ax1.set_xlabel(r'x\,[m]')
        ax1.set_ylabel(r'z\,[m]')
        ax1.set_ylim(0,z.max())
        plt.colorbar(surf,ax=axlist,format='%.0f',label=r'$p_{rel}\,[Pa]$',orientation='horizontal',pad=0.2)
        savefigPNG(animationPath + 'instantFields_' + str(i))


#%%
###############################################################################
#Creating  floor data
###############################################################################
nFloors=round(z.max()/heightFloor) #in [m]
chunkFloors=np.array_split(z, nFloors)

chunkFloorIndex=[]
for i in range(0,len(chunkFloors)):
    index=[]
    for j in range(0,len(chunkFloors[i])):
        index.append(np.argwhere(z==chunkFloors[i][j]))
    chunkFloorIndex.append(index)

chunkFloors_z0=[chunkFloors[i][0] for i in range(0,len(chunkFloors))]
chunkFloors_z1=[chunkFloors[i][-1] for i in range(0,len(chunkFloors))]
chunkFloors_index0=[chunkFloorIndex[i][0] for i in range(0,len(chunkFloorIndex))]
chunkFloors_index1=[chunkFloorIndex[i][-1] for i in range(0,len(chunkFloorIndex))]

pFBFloors=[]
pFBFloorsMean=[]
for i in range(0,len(pFB)):
    temp=[]
    temp2=[]
    for j in range(0,int(nFloors)):
        temp.append(pFB[i][int(chunkFloors_index0[j]):int(chunkFloors_index1[j]),:])
        temp2.append(np.mean(pFB[i][int(chunkFloors_index0[j]):int(chunkFloors_index1[j]),:]))
    pFBFloors.append(temp)
    pFBFloorsMean.append(temp2)

pFBFloorsMean=np.transpose(np.array(pFBFloorsMean))
FxFloorsMean=pFBFloorsMean*np.diff(z)[0]*np.diff(xFB)[0]*rho_inf

zFloors=np.arange(0,nFloors,1)


#%%Plot fx over floors

zFloors,FxFloor=createFloorData(heightFloor,z,pFB)

FxFloorsMax=np.array([np.max(FxFloor[i,:]) for i in range(0,FxFloor.shape[0])])
FxFloorsMin=np.array([np.min(FxFloor[i,:]) for i in range(0,FxFloor.shape[0])])
FxFloorsMean=np.array([np.mean(FxFloor[i,:]) for i in range(0,FxFloor.shape[0])])

latexify()
fig1, (ax1) =  plt.subplots(1)
axlist=[ax1]   
ax1.plot(FxFloorsMax/1000,zFloors,label=r'$F_{x_{max}}$')
ax1.plot(FxFloorsMean/1000,zFloors,label=r'$F_{x_{mean}}$')
ax1.plot(FxFloorsMin/1000,zFloors,label=r'$F_{x_{min}}$')
ax1.fill_betweenx(zFloors, FxFloorsMin/1000, FxFloorsMax/1000, facecolor='grey',alpha=0.2,interpolate=True)
ax1.set_xlabel(r'$F_{x_i}\,[kN]$')
ax1.set_ylabel(r'$Floors\,[-]$')
plt.legend()
savefig(rootPath + 'fx_z_floors')

#%% #%%Plot fy over floors

zFloors,FyFloor=createFloorData(heightFloor,z,pSS)
FyFloorsMax=np.array([np.max(FyFloor[i,:]) for i in range(0,FyFloor.shape[0])])
FyFloorsMin=np.array([np.min(FyFloor[i,:]) for i in range(0,FyFloor.shape[0])])
FyFloorsMean=np.array([np.mean(FyFloor[i,:]) for i in range(0,FyFloor.shape[0])])

latexify()
fig1, (ax1) =  plt.subplots(1)
axlist=[ax1]   
ax1.plot(FyFloorsMax/1000,zFloors,label=r'$F_{y_{max}}$')
ax1.plot(FyFloorsMean/1000,zFloors,label=r'$F_{y_{mean}}$')
ax1.plot(FyFloorsMin/1000,zFloors,label=r'$F_{y_{min}}$')
ax1.fill_betweenx(zFloors, FyFloorsMin/1000, FyFloorsMax/1000, facecolor='grey',alpha=0.2,interpolate=True)
ax1.set_xlabel(r'$F_{y_i}\,[kN]$')
ax1.set_ylabel(r'$Floors\,[-]$')
plt.legend()
savefig(rootPath + 'fy_z_floors')


#%%Calculate Spectra Fx
#fftX=[]
#for i in range(0,nFloors):
#    dT=np.diff(time)[0]
#    f,S= calculateSpectra(FxFloor[i,:], dT)
#    fftX.append(np.transpose(np.array([f,S])))
#
#
#
# 
#f=np.array([fftX[i][:,0] for i in range(0,len(fftX))]).T
#S=np.array([fftX[i][:,1] for i in range(0,len(fftX))]).T
#
#maxIndex=np.array(np.where(S==np.max(S)))
#fmax=f[maxIndex[0],maxIndex[1]]
#Smax=S[maxIndex[0],maxIndex[1]]
#
#latexify()
#fig1, (ax1) =  plt.subplots(1)
#axlist=[ax1]
#for i in range(0,nFloors):   
#    ax1.plot(f[:,i],S[:,i],color='grey',alpha=0.25)
#
#ax1.plot(f[:,maxIndex[1]],S[:,maxIndex[1]],color=new_colors[0],label=r'$max\{f,S,floor\}=\{'+str(fmax)+','+str(Smax)+','+str(maxIndex[1])+'\}$')
#ax1.set_xlabel(r'$f_{x_i}\,[Hz]$')
#ax1.set_ylabel(r'$S\,[-]$')
#ax1.set_xscale('log')
#ax1.set_yscale('log')
#plt.legend()
#savefig('fftx_Ampl')




#%%Calculate Spectra Fy
fftY=[]
for i in range(0,nFloors):
    dT=np.diff(time)[0]
    f,S= calculateSpectra(FyFloor[i,:], dT)
    fftY.append(np.transpose(np.array([f,S])))


f=np.array([fftY[i][:,0] for i in range(0,len(fftY))]).T
S=np.array([fftY[i][:,1] for i in range(0,len(fftY))]).T

maxIndex=np.array(np.where(S==np.max(S)))
fmax=f[maxIndex[0],maxIndex[1]]
Smax=S[maxIndex[0],maxIndex[1]]

latexify()
fig1, (ax1) =  plt.subplots(1)
axlist=[ax1]
for i in range(0,nFloors):   
    ax1.plot(f[:,i],S[:,i],color='grey',alpha=0.25)

ax1.plot(f[:,maxIndex[1]],S[:,maxIndex[1]],color=new_colors[1],label=r'$max\{f,S,floor\}=\{'+str(fmax)+','+str(Smax)+','+str(maxIndex[1])+'\}$')
ax1.set_xlabel(r'$f_{y_i}\,[Hz]$')
ax1.set_ylabel(r'$S\,[-]$')
ax1.set_xscale('log')
ax1.set_yscale('log')
plt.legend()
savefig(rootPath + 'ffty_Ampl')
