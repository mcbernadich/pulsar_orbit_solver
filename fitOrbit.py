#!/usr/bin/env python3

"""
fitOrbit.py

Script to fit orbital parameters from barycentered MJD-Ps data.
Usage:
- Give it barycentered MJD-PS data alon to plot
- Give it barycentered MJD-PS and an intial parameter file
- For details: "python3 fitOrbit.py -h"
Requires orbital_tools.py

Author: Miquel Colom i Bernadich
Date: 2026-03-16
"""

import argparse
import orbital_tools as ot

#Writen by Miquel Colom Bernadich i la mare que el va parir. Last update: 16/03/2025

parser=argparse.ArgumentParser(description="Take in barycentered MJD-Ps measurements and plot them, with or without an orbital model.")
parser.add_argument("-d","--data",help="Files with MJD, Ps and error meaurements (ms) separated by spaces. Files separated by comas, or path.")
parser.add_argument("-p","--parameter",help="Paramater file with folding ephemeris. Give only if you want to fit.")
parser.add_argument("--freq",type=bool,help="Give if the data is given in Hz instead on ms.")
args = parser.parse_args()

data=ot.dataObject(args.data,freq=args.freq)

if not args.parameter:

    print("")
    print("Will only plot the data")
    print("")
    data.plot_data()

if args.parameter:

    print("")
    params=ot.parObject(args.parameter)

    analysis=ot.analysisObject(data,params)  #It folds when it does that
    analysis.plot_folded()

    ecc_flag=input("Do you want to include eccentricity (y/n)? ")
    if ecc_flag=="n":
        print("Eccentricity will be set at 0, OM at 0, and T0=Ta.")
        params.ecc=0
        params.om=0
    fit=input("Are you ready to fit (y/n)? ")

    while fit=="n":

        par_change=input("Which paramerer would you like to modify (ps,pb,x,e,om,t0)? ")
        if par_change=="ps":
            params.f0=1000/float(input("Input new Ps (ms): "))
        if par_change=="pb":
            params.pb=float(input("Input new Pb (days): "))
        if par_change=="x":
            params.x=float(input("Input new x (ls): "))
        if par_change=="e":
            params.ecc=float(input("Input new ecc: "))
        if par_change=="om":
            params.om=float(input("Input new OM (deg): "))
        if par_change=="t0":
            params.t0=float(input("Input new T0 (TA if circular) (days): "))
        if ecc_flag=="n":
            params.ecc=0
            params.om=0

        print("")
        print("- Pulsar frequency: {} Hz".format(params.f0))
        print("- Pulsar spin: {} ms".format(1000/params.f0))
        print("- Orbital period: {} days".format(params.pb))
        print("- Projected axis: {} ls".format(params.x))
        print("- Eccentricity: {}".format(params.ecc))
        print("- Periastron angle: {} degrees".format(params.om))
        print("- Periastron passage: {} (MJD)".format(params.t0))
        print("")

        analysis.refold()
        analysis.plot_folded()


        fit=input("Are you ready to fit (y/n)?")
        print("")
    
    if ecc_flag=="n":

        analysis.fit_circular()

    elif ecc_flag=="y":

        analysis.fit_eccentric()

    analysis.refold()
    analysis.plot_publication()
    
    new_file=input("Give a file name to save the solution: ")
    params.write(new_file)