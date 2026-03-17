#!/usr/bin/env python3

"""
searchOrbit.py

Script to search for the best orbital period to a set of barycentered MJD-Ps data.
Usage:
- Give it barycentered MJD-PS and an intial parameter
- Give it an orbital period ragne to search and, optionally, an eccentricity range
- Intended to find the best orbital period, it does not report a full Keplerian fit. Use fitOrbit.py once you have the orbital period.
- For details: "python3 searchOrbit.py -h"
Requires orbital_tools.py

Author: Miquel Colom i Bernadich
Date: 2026-03-16
"""

import argparse
import orbital_tools as ot
import numpy as np
import time
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore", category=SyntaxWarning)

parser=argparse.ArgumentParser(description="Searcht the correct orbital from MJD-Ps measurements.")
parser.add_argument("-d","--data",help="Files with MJD, Ps and error meaurements (ms) separated by spaces. Files separated by comas, or path.")
parser.add_argument("-p","--parameter",help="Paramater file with folding ephemeris. Orbital period and eccentricity are ignored, but not the other parameters (you just need reasonable guesses).")
parser.add_argument("-e","--eccentricity",help="Eccentricity range to search for, {min}:{max}. Fit for eccentricity and periastron. Default: circular orbit")
parser.add_argument("-r","--range",help="Period range to search, {min}:{max}, in days.")
parser.add_argument("-x","--chi2r",type=float,default=10,help="Maximum reduced chi2 for solutions (default: <10).")
parser.add_argument("--freq",type=bool,help="Give if the data is given in Hz instead on ms.")
args = parser.parse_args()

if args.eccentricity:
    ecc_flag="y"
    ecc_lower=float(args.eccentricity.split(":")[0])
    ecc_upper=float(args.eccentricity.split(":")[1])
else:
    ecc_flag="n"

data=ot.dataObject(args.data,freq=args.freq)
print("")
params=ot.parObject(args.parameter)
analysis=ot.analysisObject(data,params)

trial_pb=float(args.range.split(":")[0])
last_pb=float(args.range.split(":")[1])
data_span=data.time[-1]-data.time[0]
print("The data span is",data_span,"days.")

chi2rs=[]
trial_pbs=[]
peak_pbs=[]
peak_chi2rs=[]
i=0
correction_attempts=0
unsolvable_knots=0

start=time.time()

# Loop of Pb trials in the specified bounds
while trial_pb <= last_pb:

    if np.mod(i,1000)==0 and i==0:
        print("Current trial Pb: {} days. Current step: {} days.".format(trial_pb,1e-2*(trial_pb**2)/(2*np.pi*data_span)))
    if np.mod(i,1000)==0 and i!=0:
        print("Current trial Pb: {} days. Current step: {} days. Lowest CHI2R so far: {}".format(trial_pb,1e-2*(trial_pb**2)/(2*np.pi*data_span),round(np.min(chi2rs),3)))

    params.pb=trial_pb
    analysis.refold()

    step=1e-2*(trial_pb**2)/(data_span) #The last point only moves by an orbital phase of 0.01

    #Fit circular orbits if eccentricity not specified
    if ecc_flag=="n":
    
        bounds=([0,trial_pb-step*0.5,0,-np.inf],[np.inf,trial_pb+step*0.5,np.inf,np.inf])

        try:
            analysis.fit_circular(print_result=False,bounds=bounds)
            chi2rs.append(analysis.chi2r)
        except:
            print("Failed at Pb=",trial_pb,"days, i=",i)
            chi2rs.append(np.max(chi2rs))
        trial_pbs.append(trial_pb)

    #Fit eccentric orbits
    elif ecc_flag=="y":

        bounds=([0,trial_pb-step*0.5,0,ecc_lower,0,-np.inf],[np.inf,trial_pb+step*0.5,np.inf,ecc_upper,360,np.inf])

        try:
            analysis.fit_eccentric(print_result=False,bounds=bounds)
            chi2rs.append(analysis.chi2r)
        except:
            print("Failed at Pb=",trial_pb,"days, i=",i)
            chi2rs.append(np.max(chi2rs))
        trial_pbs.append(trial_pb)

    # Store minima
    if i>1 and chi2rs[-1]>chi2rs[-2] and chi2rs[-2]<chi2rs[-3] and chi2rs[-2]<args.chi2r:
        print("Peak found at {} days, with CHI2R={}".format(trial_pbs[-2],round(chi2rs[-2],3)))
        peak_chi2rs.append(chi2rs[-2])
        peak_pbs.append(trial_pbs[-2])

    trial_pb=trial_pb+step
    i=i+1

print("Run time: {} s".format(time.time()-start))

chi2rs=np.array(chi2rs)
trial_pbs=np.array(trial_pbs)
print("Total number of trials: {}".format(np.size(chi2rs)))
print("Number of fits with chi2r<{}: {}".format(args.chi2r,np.size(peak_chi2rs)))

best_pb=trial_pbs[np.argmin(chi2rs)]

print("")
print("Best orbital period= {} days".format(best_pb))
print("Best chi2r= {}".format(np.min(chi2rs)))
print("")

plt.plot(trial_pbs,chi2rs)
plt.plot([],[]," ",label=r"$P_\mathrm{b}=$ "+str(round(best_pb,5))+" d")
plt.ylabel(r"$\chi^2_\mathrm{r}$",size=14)
plt.xlabel("Trial orbital period (days)",size=14)
plt.xlim(float(args.range.split(":")[0]),float(args.range.split(":")[1]))
plt.tick_params(labelsize=12)
plt.title(params.psr_name,size=14)
plt.yscale("log")
plt.legend(fontsize=12)
plt.tight_layout()
#plt.savefig("chi2rs.png",format="png",dpi=800)
plt.show()

np.savetxt(params.psr_name+".out",np.array([trial_pbs,chi2rs]).T,header="trial Pb (days), reduced chi2\n------------------------------")
np.savetxt(params.psr_name+".best",np.array([best_pb,np.min(chi2rs)]).T,header="trial Pb (days), reduced chi2\n------------------------------")
np.savetxt(params.psr_name+".peaks",np.array([peak_pbs,peak_chi2rs]).T,header="trial Pb (days), reduced chi2\n------------------------------")