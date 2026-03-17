#!/usr/bin/env python3

"""
searchOrbit.py

Script to estimate orbital parameters from Ps-acc measurements (pdmp), Ps-Pdot (prepfols) or f0-f1 (tempo/tempo2)
Usage:
- Give it barycentered Ps-acc data 
- Fit for a circular or eccentric orbit 
- well-tested with moderate eccentricity, can't promise for e~0.9
- For details: "python3 fitEllise.py -h"
Requires orbital_tools.py

Author: Miquel Colom i Bernadich
Date: 2026-03-17
"""

import numpy as np
import orbital_tools as ot
import argparse
import warnings
warnings.filterwarnings("ignore", category=SyntaxWarning)
warnings.filterwarnings("ignoren", category=RuntimeWarning)

parser=argparse.ArgumentParser(description="Take in Ps-Psdot measurements and fit orbital parameters.")
parser.add_argument("-d","--data",help="File with MJD, Ps (ms), accelerations and error (m/s2) separated by spaces (pdmp mode).")
parser.add_argument("-f","--freq",type=bool,help="Set if data is given in Hz and Hz/s instead on ms and ms/s2 (tempo/tempo2 mode).")
parser.add_argument("-s","--spin",type=bool,help="Set if data given in ms and s/s instead on ms and ms/2 (prepfold mode).")
parser.add_argument("-n","--name",help="Name of the pulsar")
args = parser.parse_args()

ellipse_data=ot.ellipseObject(args.data,freq=args.freq)
ellipse_data.fit_ellipse_coarse()
ellipse_data.plot_ellipse_circular(args.name)

print("")
print("Ready to make a more serious fit.")
if np.size(ellipse_data.period)>=5:
    ecc_flag=input("There are at least 5 data points. Do you want to include eccentricity (y/n)? ")
else:
    ecc_flag="n"

if ecc_flag=="n":

    ellipse_data.fit_ellipse_circular()
    ellipse_data.plot_ellipse_circular(args.name)

if ecc_flag=="y":

    ellipse_data.fit_ellipse_eccentric()
    ellipse_data.plot_ellipse_eccentric(args.name)
