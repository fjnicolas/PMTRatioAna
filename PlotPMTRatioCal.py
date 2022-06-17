import numpy as np
from matplotlib import pyplot as plt
import math
import ROOT
plt.rcParams.update(plt.rcParamsDefault)
params = {'legend.fontsize': 'xx-large',
          'figure.figsize': (10, 7),
         'axes.labelsize': 'xx-large',
         'axes.titlesize':'xx-large',
         'xtick.labelsize':'large',
         'ytick.labelsize':'large'}
plt.rcParams.update(params)
import pylab
import os, glob, sys
import seaborn as sns
import pandas as pd
import json
from matplotlib import gridspec
from matplotlib.ticker import (AutoMinorLocator, MultipleLocator)
import argparse
import seaborn as sns
import matplotlib.lines as mlines
import uproot



parser = argparse.ArgumentParser()
parser.add_argument("-s", "--Filepath", help="Input file name",  default="none")
parser.add_argument("-n", "--MaxEvents", help="Max Number of Events",  default=-1, type=int)
parser.add_argument("-nfiles", "--NFiles", help="NFiles",  default=1e4, type=int)
parser.add_argument("-cut", "--CutOption", help="Option cut",  default=0, type=int)
parser.add_argument("-origin", "--Origin", help="Origin",  default=0, type=int)
parser.add_argument("-sat", "--Saturation", help="Exclude saturated channels",  default=0, type=int)
parserargs = parser.parse_args()

InputFileList=[]
for filepath in glob.iglob(parserargs.Filepath):
    InputFileList.append(filepath)



MyCut=""
if(parserargs.CutOption==0):
    print("No cuts...", MyCut)

MyCut=MyCut+"TotalEnergyDep>=100"
print("CUT:", MyCut)


pdtype=[];
f=open(os.getenv('WD')+"Utilities/GenericDisplays/include/sbndPDSMap_v02_00.txt","r");
lines=f.readlines();
f.close()
for x in lines:
    _pdname=x.split(' ')[4]
    if '\n' in _pdname: _pdname=_pdname[:-1]
    pdtype.append(_pdname)



def ReadVars(fTree, varnames="", cut=""):
    cutF=ROOT.TCut(cut)
    h2 = fTree.Draw(varnames, MyCut)
    X=list(  np.ndarray((h2), 'd', fTree.GetV2()) )
    PMTRatio=list(np.ndarray((h2), 'd', fTree.GetV1()))
    return PMTRatio, X




XSimPhotons=[]
PMTRatioSimPhotons=[]
PMTRatioPE=[]
XPE=[]

varnamesMC="PMTRatioSimPhotons:dEPromX"
varnamesReco="PMTRatioPE:dEPromX"


for file_ix, filepath in enumerate(InputFileList):
    if(file_ix>parserargs.NFiles): continue
    print("READING FILE: ", file_ix, filepath)

    f = ROOT.TFile.Open(filepath)
    tree =  f.Get("pmtratioana/PMTRatioTree")
    print(tree)

    v1, v2 = ReadVars(tree, varnamesMC, MyCut)
    PMTRatioSimPhotons.extend(v1)
    XSimPhotons.extend(v2)

    v1, v2 = ReadVars(tree, varnamesReco, MyCut)
    PMTRatioPE.extend(v1)
    XPE.extend(v2)


fig = plt.figure(figsize=(14, 8), num=str(MyCut))
params = {'legend.fontsize': 20,
          'figure.figsize': (10, 7),
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
plt.rcParams.update(params)
ax1=fig.add_subplot(221); ax1.grid()
ax2=fig.add_subplot(222); ax2.grid()
ax3=fig.add_subplot(223); ax3.grid()
ax4=fig.add_subplot(224); ax4.grid()

XSimPhotons=np.abs(np.array(XSimPhotons))
PMTRatioSimPhotons=np.array(PMTRatioSimPhotons)
ax1.scatter(PMTRatioSimPhotons, XSimPhotons, label="MC")
ax1.set_xlabel("PMTRatio"); ax1.set_ylabel(r"$<X>_{MC}$ [cm]")
ax1.legend()

XPE=np.abs(np.array(XSimPhotons))
PMTRatioPE=np.array(PMTRatioSimPhotons)
ax2.scatter(PMTRatioPE, XPE, label="Reco")
ax2.set_xlabel("PMTRatio"); ax2.set_ylabel(r"$<X>_{MC}$ [cm]")
ax2.legend()

plt.show()
