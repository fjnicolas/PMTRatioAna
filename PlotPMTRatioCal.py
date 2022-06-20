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
MyCut=MyCut+"&& dESpreadX<5"
MyCut=MyCut+"&& dETPCBalance>0.99"
MyCut=MyCut+"&& abs(dEPromY)<180 && dEPromZ>25 && dEPromZ<475"
print("CUT:", MyCut)

def GetPMTRatioProfile( Xv, Yv, bin_range=[0., 1.], binS=0.025):
    NBinsF=int( (bin_range[1]-bin_range[0]) / binS )
    hP=ROOT.TProfile("","",NBinsF, bin_range[0], bin_range[1], "s");

    for j in range(len(Xv)):
        hP.Fill(Xv[j], Yv[j])

    Xp=[]
    Yp=[]
    YpErr=[]
    for ix in range(1, hP.GetNbinsX()):
        Xp.append(hP.GetBinCenter(ix))
        Yp.append(hP.GetBinContent(ix))
        YpErr.append(hP.GetBinError(ix))

    return hP, [Xp, Yp, YpErr]

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


fig = plt.figure(figsize=(13, 8), num=str(MyCut))
params = {'legend.fontsize': 20,
          'figure.figsize': (10, 7),
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
plt.rcParams.update(params)
ax1=fig.add_subplot(211); ax1.grid()
ax2=fig.add_subplot(212); ax2.grid()
#ax3=fig.add_subplot(223); ax3.grid()
#ax4=fig.add_subplot(224); ax4.grid()


binFSize=0.01
minF=0.04
maxF=0.8
Alpha=0.3


XSimPhotons=np.abs(np.array(XSimPhotons))
PMTRatioSimPhotons=np.array(PMTRatioSimPhotons)
hp_SimPhotons, v_SimPhotons = GetPMTRatioProfile(PMTRatioSimPhotons, XSimPhotons, [minF, maxF], binFSize)
ax1.scatter(PMTRatioSimPhotons, XSimPhotons, marker="o", s=0.2, alpha=Alpha, c='gray')
ax1.errorbar(v_SimPhotons[0], v_SimPhotons[1], v_SimPhotons[2], fmt="v", label="MC", marker="v", color="C0", markersize=5, linewidth=2, capsize=6)
ax1.set_xlabel("PMTRatio"); ax1.set_ylabel(r"$<X>$ [cm]")
ax1.legend()

XPE=np.abs(np.array(XSimPhotons))
PMTRatioPE=np.array(PMTRatioSimPhotons)
hp_PE, v_PE = GetPMTRatioProfile(PMTRatioPE, XPE, [minF, maxF], binS=binFSize)
ax2.scatter(PMTRatioPE, XPE, marker="o", s=0.2, alpha=Alpha, c='gray')
ax2.errorbar(v_PE[0], v_PE[1], v_PE[2], fmt="^", label="Reco", marker="^", color="C1", markersize=5, linewidth=2, capsize=6)
ax2.set_xlabel("PMTRatio"); ax2.set_ylabel(r"$<X>$ [cm]")
ax2.legend()


fileOut = ROOT.TFile("PMTRatioCalibration.root", "RECREATE");
fileOut.WriteObject(hp_PE, "PMTRatioCalibrationProfile");
fileOut.Close();


plt.show()
