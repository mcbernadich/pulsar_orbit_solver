#!/usr/bin/env python3

"""
orbital_tools.py

Utilities to fit for orbits.
Provides:
- Orbital model functions: Kepler equations, circular/eccentric period prediction, etc.
- dataObject: load and manage spin period measurements
- parObject: read/write pulsar ephemeris
- analysisObject: fold, fit, and plot
- ellipseObject: read data and fit the acc-Ps ellipse to estiamte orbital parameters

Author: Miquel Colom i Bernadich
Date: 2026-03-17
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import newton
from scipy.optimize import curve_fit
from scipy.optimize import least_squares
from scipy.optimize import differential_evolution

def kepler_equation(E,M,e):
    return E-e*np.sin(E)-M

def kepler_equation_derivative(E,M,e):
    return 1-e*np.cos(E)

#Return the baricentered spin frequency from the eccentric anomaly
def frequency_fromE(E,Edot,f0,fxom,fxomecc):
    return f0+fxom*np.sin(E)*Edot-fxomecc*np.cos(E)*Edot #This can be vectorial too.

#Return the barycentered frequency derivatve from the eccentric anomaly
def frequency_derivative_fromE(E,Edot,Edotdot,fxom,fxomecc):
    return fxom*(np.sin(E)*Edotdot+np.cos(E)*np.square(Edot))-fxomecc*(np.cos(E)*Edotdot-np.sin(E)*np.square(Edot)) #This one can be vectorial too.

#Return ps as a function of time for eccentric orbits
def predict_period(time,f0,pb,x,e,omega,t0):

    #Compute eccentric anomaly and time derivatives.
    M=2*np.pi*(time-t0)/pb
    E=newton(kepler_equation,M,fprime=kepler_equation_derivative,args=(M,abs(e)))
    pi2orb=2*np.pi/(pb*86400)
    Edot=pi2orb/(1-abs(e)*np.cos(E))
    Edotdot=-(abs(e)*np.sin(E)*Edot**2)/(1-abs(e)*np.cos(E))

    #Compute constants for the frequency equations.
    fxom=f0*abs(x)*np.sin(np.deg2rad(omega))
    fxomecc=f0*abs(x)*np.cos(np.deg2rad(omega))*np.sqrt(1-e**2)

    #Compute periods.
    ps=1/frequency_fromE(E,Edot,f0,fxom,fxomecc)

    return ps

#Return ps as a function of time for circular orbits
def predict_period_circular(time,f0,pb,x,t0):

    #Compute eccentric anomaly and time derivatives.
    E=2*np.pi*(time-t0)/pb
    pi2orb=2*np.pi/(pb*86400)
    Edot=pi2orb
    Edotdot=0

    #Compute constants for the frequency equations.
    fxom=f0*abs(x)*np.sin(np.deg2rad(0))
    fxomecc=f0*abs(x)*np.cos(np.deg2rad(0))

    #Compute periods.
    ps=1/frequency_fromE(E,Edot,f0,fxom,fxomecc)

    return ps

#Return the fdot given an f and orbital parameters in circular orbits
def fdot_from_f_circular(f,ps,pb,x):

    x=abs(x)
    pb=abs(pb)*86400
    cosE=(1-f*ps)*pb/2/np.pi/x
    return ((1/ps)*x*np.sqrt(1-cosE**2)*(2*np.pi/pb)**2,-(1/ps)*x*np.sqrt(1-cosE**2)*(2*np.pi/pb)**2)

#Compute the residulas (actually chi2 before sum) given data and the predictiona bove 
def residuals_circular(pars,p,pdot,errors):

    (ps,pb,x)=(pars[0],pars[1],pars[2])
    (model_pos,model_neg)=fdot_from_f_circular(1/p,ps,pb,x)

    #Transform into pdot values.
    model_pos=-model_pos*p**2
    model_neg=-model_neg*p**2

    residuals=np.minimum(np.abs(model_pos-pdot),np.abs(model_neg-pdot))
    residuals[residuals==np.nan]=999999999999    

    return (residuals/errors)**2

#Return the fdot given an f and orbital parameters in eccentric orbits
def fdot_from_f_eccentric(f,ps,pb,x,e,om):

    x=abs(x)
    e=abs(e)
    om=om%360
    pb=abs(pb)*86400
    pi2orb=2*np.pi/pb
    fr=f-1/ps
    fxom=pi2orb*abs(x)*np.sin(np.deg2rad(om))/ps
    fxomecc=pi2orb*abs(x)*np.cos(np.deg2rad(om))*np.sqrt(1-e**2)/ps

    #It is no longer symeetric due to eccentricity. E here represents the eccentric anomaly at positive and negative fdots
    cosE1=(fr*fxomecc+e*fr**2-fxom*np.sqrt(fxom**2+fxomecc**2+2*e*fxomecc*fr+(e**2-1)*fr**2))/(fxom**2+fxomecc**2+2*e*fxomecc*fr+e**2*fr**2)
    cosE2=(fr*fxomecc+e*fr**2+fxom*np.sqrt(fxom**2+fxomecc**2+2*e*fxomecc*fr+(e**2-1)*fr**2))/(fxom**2+fxomecc**2+2*e*fxomecc*fr+e**2*fr**2)

    fr=fr+1e-12 #Avoid numerical issues at fr=0

    dcosE1=(fr*fxomecc+e*fr**2-fxom*np.sqrt(fxom**2+fxomecc**2+2*e*fxomecc*fr+(e**2-1)*fr**2))/(fxom**2+fxomecc**2+2*e*fxomecc*fr+e**2*fr**2)-cosE1
    dcosE2=(fr*fxomecc+e*fr**2+fxom*np.sqrt(fxom**2+fxomecc**2+2*e*fxomecc*fr+(e**2-1)*fr**2))/(fxom**2+fxomecc**2+2*e*fxomecc*fr+e**2*fr**2)-cosE2

    sinE1=-np.sqrt(1-cosE1**2)
    sinE1[dcosE1<0]=np.sqrt(1-cosE1[dcosE1<0]**2)
    sinE2=np.sqrt(1-cosE2**2)
    sinE2[dcosE2<0]=-np.sqrt(1-cosE2[dcosE2<0]**2)

    Edot1=pi2orb/(1-e*cosE1)
    Edot2=pi2orb/(1-e*cosE2)
    Edotdot1=-(e*sinE1*Edot1**2)/(1-e*cosE1)
    Edotdot2=-(e*sinE2*Edot2**2)/(1-e*cosE2)

    fdot1=fxom*(sinE1*Edotdot1+cosE1*np.square(Edot1))/pi2orb-fxomecc*(cosE1*Edotdot1-sinE1*np.square(Edot1))/pi2orb
    fdot2=fxom*(sinE2*Edotdot2+cosE2*np.square(Edot2))/pi2orb-fxomecc*(cosE2*Edotdot2-sinE2*np.square(Edot2))/pi2orb

    return fdot1,fdot2

#Compute the residulas (actually chi2 before sum) given data and the predictiona bove 
def residuals_eccentric(pars,p,pdot,errors):

    (ps,pb,x,e,om)=(pars[0],pars[1],pars[2],pars[3],pars[4])
    (model1,model2)=fdot_from_f_eccentric(1/p,ps,pb,x,e,om)

    #Transform into pdot values.
    model1=-model1*p**2
    model2=-model2*p**2

    residuals=np.minimum(np.abs(model1-pdot),np.abs(model2-pdot))
    residuals[np.isnan(residuals)]=999999999999

    return np.sum((residuals/errors)**2)

# ----------------------
# dataObject: time-series
# ----------------------
class dataObject:
    def __init__(self, path, freq=False):

        import glob

        # Expand wildcards or split comma-separated list
        if "*" in path or "?" in path:
            files = sorted(glob.glob(path))
        else:
            files = [f.strip() for f in path.split(",")]

        if not files:
            raise ValueError(f"No files found matching {path}")

        # Load and concatenate
        arrays = [np.loadtxt(f) for f in files]
        self.data = np.concatenate(arrays).T

        # Optional frequency transform
        if freq:
            self.apply_frequency_transform()

        # Convenience attributes
        self.time = self.data[0]
        self.period = self.data[1]
        self.error = self.data[2]

    def apply_frequency_transform(self):
        #Convert Hz to ms
        self.data[2] = 1000 * self.data[2] / (self.data[1] ** 2)
        self.data[1] = 1000 / self.data[1]

    def plot_data(self,show=True):
        plt.ylabel(r"$P_\mathrm{SSB}$ (ms)",size=12)
        plt.xlabel("Time (MJD)",size=12)
        plt.errorbar(self.time,self.period,yerr=self.error,fmt="bo",label="data",markersize=5)
        plt.grid()
        plt.legend()
        plt.tight_layout()

        if show:
        	plt.show()

# ----------------------
# parObject: ephemeris
# ----------------------
class parObject:
    def __init__(self, file):
        """
        Read spin and Keplerian parameter from file.
        """
        # Initialize all expected attributes
        self.file = file
        self.psr_name = None
        self.f0 = None
        self.df0 = 0
        self.pb = None
        self.dpb = 0
        self.x = None
        self.dx = 0
        self.ecc = 0
        self.decc = 0 
        self.om = 0
        self.dom = 0
        self.t0 = None
        self.dt0 = 0

        # Read ephemeris file
        with open(self.file, "r") as ephemeris:
            for line in ephemeris:
                line = line.strip().split()
                if not line:
                    continue  # skip empty lines
                key = line[0].upper()

                if key == "PSR" or key == "PSRJ":
                    self.psr_name = line[1]
                if key == "F0":
                    self.f0 = float(line[1])
                    try:
                        self.df0 = float(line[3])
                    except (IndexError, ValueError):
                        self.df0 = 0
                elif key == "P0":
                    self.f0 = 1 / float(line[1])  # convert period to frequency
                    try:
                        self.df0 = float(line[3]) / float(line[1])**2
                    except (IndexError, ValueError):
                        self.df0 = 0
                elif key == "PB":
                    self.pb = float(line[1])
                    try:
                        self.dpb = float(line[3])
                    except (IndexError, ValueError):
                        self.dpb = 0
                elif key == "A1":
                    self.x = float(line[1])
                    try:
                        self.dx = float(line[3])
                    except (IndexError, ValueError):
                        self.dx = 0
                elif key == "ECC":
                    self.ecc = float(line[1])
                    try:
                        self.decc = float(line[3])
                    except (IndexError, ValueError):
                        self.decc = 0
                elif key == "OM":
                    self.om = float(line[1])
                    try:
                        self.dom = float(line[3])
                    except (IndexError, ValueError):
                        self.dom = 0
                elif key == "T0":
                    self.t0 = float(line[1])
                    try:
                        self.dt0 = float(line[3])
                    except (IndexError, ValueError):
                        self.dt0 = 0

        # Check that all required parameters are set
        if None in (self.f0, self.pb, self.x, self.t0):
            sys.exit("The ephemeris file doesn't include all the necessary parameters.")
        
        print("Pulsar parameters loaded from {}:".format(file))
        print("- Pulsar frequency: {} Hz".format(self.f0))
        print("- Pulsar spin: {} ms".format(1000/self.f0))
        print("- Orbital period: {} days".format(self.pb))
        print("- Projected axis: {} ls".format(self.x))
        print("- Eccentricity: {}".format(self.ecc))
        print("- Periastron angle: {} degrees".format(self.om))
        print("- Periastron passage: {} (MJD)".format(self.t0))
        print("")

    def write(self, new_file):
        """
        Write the current ephemeris values to a new file.
        """
        with open(self.file, "r") as ephemeris, open(new_file, "w") as out:
            for line_str in ephemeris:
                line = line_str.strip().split()
                if not line:
                    out.write(line_str)
                    continue

                key = line[0].upper()  # "F0" -> "F0"

                if key in ("F0", "P0"):
                    # write P0 using f0
                    out.write(f"P0 {1/self.f0} 1 {self.df0/self.f0**2}\n")
                elif key == "PB":
                    out.write(f"PB {self.pb} 1 {self.dpb}\n")
                elif key == "A1":
                    out.write(f"A1 {self.x} 1 {self.dx}\n")
                elif key == "ECC":
                    out.write(f"ECC {self.ecc} 1 {self.decc}\n")
                elif key == "OM":
                    out.write(f"OM {self.om} 1 {self.dom}\n")
                elif key == "T0":
                    out.write(f"T0 {self.t0} 1 {self.dt0}\n")
                else:
                    out.write(line_str)

        return new_file

# ----------------------
# analysisObject: of analysis
# ----------------------
class analysisObject:
    def __init__(self,data,par):
        self.data = data
        self.par = par

        self.folded_time = None
        self.folded_t0 = None

        self.fold()

    def fold(self):
        """
        Reduce the time-series data modulo the orbital period of the ParObject.
        Returns a new array with folded times and associated period/error.
        """

        par = self.par

        if par.pb is None:
            raise ValueError("ParObject must have pb defined for reduction")

        time = self.data.time

        self.folded_time = ( time - time[0] ) % par.pb
        self.folded_t0 = ( par.t0 - time[0] ) % par.pb

    def refold(self):
        """
        Recompute quantities that depend on par.
        """
        self.fold()

    def plot_folded(self,show=True):

        par = self.par
        data = self.data

        times=np.arange(-par.pb/20,par.pb+par.pb/20,par.pb/500)
        modelps=predict_period(times,par.f0,par.pb,par.x,par.ecc,par.om,self.folded_t0)
        modelps=1000*modelps        

        plt.ylabel(r"$P_\mathrm{SSB}$ (ms)",size=12)
        plt.xlabel("Time (MJD)",size=12)
        plt.plot(times,modelps,"c-",label="model")
        plt.errorbar(self.folded_time,data.period,yerr=data.error,fmt="bo",label="data",markersize=5)
        plt.vlines(self.folded_t0,np.min(modelps)-(np.max(modelps)-np.min(modelps))/20,np.max(modelps)+(np.max(modelps)-np.min(modelps))/20,colors="r",linestyles="dashed",label="periastron")
        plt.hlines(1000/par.f0,-par.pb/20,par.pb+par.pb/20,colors="g",linestyles="dotted",label="true spin")
        plt.xlim(-par.pb/20,par.pb+par.pb/20)
        plt.ylim(np.min(modelps)-(np.max(modelps)-np.min(modelps))/20,np.max(modelps)+(np.max(modelps)-np.min(modelps))/20)
        plt.grid()
        plt.legend()
        plt.tight_layout()

        if show:
            plt.show()

    def plot_publication(self,show=True):

        par = self.par
        data = self.data

        times=np.arange(-par.pb/20,par.pb+par.pb/20,par.pb/500)
        modelps=predict_period(times,par.f0,par.pb,par.x,par.ecc,par.om,self.folded_t0)
        modelps=1000*modelps        

        plt.ylabel(r"$P_\mathrm{SSB}$ (ms)",size=14)
        plt.xlabel("Time (MJD)",size=14)
        plt.plot(times/par.pb,modelps,"c-",label="model")
        plt.errorbar(self.folded_time/par.pb,data.period,yerr=data.error,fmt="bo",label="data",markersize=5)
        plt.xlim(-1/20,1+1/20)
        plt.ylim(np.min(modelps)-(np.max(modelps)-np.min(modelps))/20,np.max(modelps)+(np.max(modelps)-np.min(modelps))/20)
        plt.title(f"PSR {par.psr_name}",size=14)
        plt.tick_params(labelsize=12)
        plt.tight_layout()

        if show:
            plt.show()

    def fit_circular(self,print_result=True,bounds=None):

        par = self.par
        data = self.data

        if bounds:
            (fit_par, cov) = curve_fit(predict_period_circular, data.time, data.period/1000, p0=[par.f0,par.pb,par.x,par.t0], sigma=data.error/1000, absolute_sigma=True,bounds=bounds)
        else:
            (fit_par, cov) = curve_fit(predict_period_circular, data.time, data.period/1000, p0=[par.f0,par.pb,par.x,par.t0], sigma=data.error/1000, absolute_sigma=True)
        (par.f0,par.pb,par.x,par.ecc,par.om,par.t0) = (fit_par[0],fit_par[1],abs(fit_par[2]),0,0,fit_par[3])
        (par.df0,par.dpb,par.dx,par.decc,par.dom,par.dt0) = (np.sqrt(cov[0,0]),np.sqrt(cov[1,1]),np.sqrt(cov[2,2]),0,0,np.sqrt(cov[3,3]))
        chi2=np.sum(((data.period-1000*predict_period_circular(data.time,par.f0,par.pb,par.x,par.t0))/data.error)**2)
        self.chi2r=chi2/(np.size(data.time)-4)
        
        if print_result==True:
            print("- Pulsar frequency: {}+-{} Hz".format(par.f0,par.df0))
            print("- Pulsar spin: {}+-{} ms".format(1000/par.f0,1000*par.df0/(par.f0**2)))
            print("- Orbital period: {}+-{} days".format(par.pb,par.dpb))
            print("- Projected axis: {}+-{} ls".format(abs(par.x),par.dx))
            print("- Time of ascending node: {}+-{} (MJD)".format(par.t0,par.dt0))
            print("- CHI2R: {}".format(self.chi2r))
            print("")

    def fit_eccentric(self,print_result=True,bounds=None):

        par = self.par
        data = self.data

        if bounds:
            (fit_par, cov) = curve_fit(predict_period, data.time, data.period/1000, p0=[par.f0,par.pb,par.x,par.ecc,par.om,par.t0], sigma=data.error/1000, absolute_sigma=True,bounds=bounds)
        else:
            (fit_par, cov) = curve_fit(predict_period, data.time, data.period/1000, p0=[par.f0,par.pb,par.x,par.ecc,par.om,par.t0], sigma=data.error/1000, absolute_sigma=True)
        (par.f0,par.pb,par.x,par.ecc,par.om,par.t0) = (fit_par[0],fit_par[1],abs(fit_par[2]),abs(fit_par[3]),np.mod(fit_par[4],360),fit_par[5])
        (par.df0,par.dpb,par.dx,par.decc,par.dom,par.dt0) = (np.sqrt(cov[0,0]),np.sqrt(cov[1,1]),np.sqrt(cov[2,2]),np.sqrt(cov[3,3]),np.sqrt(cov[4,4]),np.sqrt(cov[5,5]))
        chi2=np.sum(((data.period-1000*predict_period(data.time,par.f0,par.pb,par.x,par.ecc,par.om,par.t0))/data.error)**2)
        self.chi2r=chi2/(np.size(data.time)-6)

        if print_result==True:
            print("- Pulsar frequency: {}+-{} Hz".format(par.f0,par.df0))
            print("- Pulsar spin: {}+-{} ms".format(1000/par.f0,1000*par.df0/(par.f0**2)))
            print("- Orbital period: {}+-{} days".format(par.pb,par.dpb))
            print("- Projected axis: {}+-{} ls".format(abs(par.x),par.dx))
            print("- Eccentricity: {}+-{}".format(abs(par.ecc),par.decc))
            print("- Periastron angle: {}+-{} degrees".format(np.mod(par.om,360),par.dom))
            print("- Time of periastron: {}+-{} (MJD)".format(par.t0,par.dt0))
            print("- CHI2R: {}".format(self.chi2r))
            print("")

# ----------------------
# ellipseObject: load and analyse ellipse data
# ----------------------
class ellipseObject:
    def __init__(self, path, freq=False, spin=False):

        import glob

        # Expand wildcards or split comma-separated list
        if "*" in path or "?" in path:
            files = sorted(glob.glob(path))
        else:
            files = [f.strip() for f in path.split(",")]

        if not files:
            raise ValueError(f"No files found matching {path}")

        # Load and concatenate
        arrays = [np.loadtxt(f) for f in files]
        self.data = np.concatenate(arrays).T

        # Optional frequency transform
        if freq:
            self.apply_frequency_transform()
        elif spin:          # Optionalspin transform
            self.apply_spin_transform()
        else:               # Convert ms to s
            self.period=self.data[1]/1000
            self.acceleration=self.data[2]
            self.acc_error=self.data[3]
            self.derivative=self.period*self.acceleration/299792458
            self.der_error=self.period*self.acc_error/299792458

        # Convenience attributes (unused)
        self.time = self.data[0]

    def apply_frequency_transform(self):
        #Convert Hz anf Hz/s to s and s/s
        self.der_error=self.data[3]/(self.data[1]**2)
        self.derivative=-self.data[2]/(self.data[1]**2)
        self.period=1/self.data[1]
        #Convert s/s to m/s2
        self.acceleration=299792458*self.derivative/self.period
        self.acc_error=299792458*self.der_error/self.period

    def apply_spin_transform(self):
        #Convert ms and s/s s and m/s2
        self.period=self.data[1]/1000
        self.derivative=self.data[2]
        self.der_error=self.data[3]
        self.acceleration=299792458*self.derivative/self.period
        self.acc_error=299792458*self.der_error/self.period

    def fit_ellipse_coarse(self):

        print("")
        print("Assuming circular and making a first, coarse fit...")
        print("")

        #Do a first crude fit without uncertainties.
        
        coeff=np.polynomial.polynomial.polyfit(self.period,self.acceleration**2,2)
        (a0,a1,a2)=(coeff[2],coeff[1],coeff[0])

        self.ps=-a1/(2*a0)
        self.pb=2*np.pi*299792458/(self.ps*np.sqrt(-a0))
        self.x=self.pb*np.sqrt(self.ps**2-a2/a0)/(2*np.pi*self.ps)
        self.acc_amp=np.square(2*np.pi/self.pb)*self.x*299792458
        self.ps_amp=2*np.pi*self.ps*self.x/self.pb
        self.pb=self.pb/24/3600

        residuals=residuals_circular([self.ps,self.pb,self.x],self.period,self.derivative,self.der_error)

        chi2=np.sum(residuals**2)
        self.chi2r=chi2/(np.size(self.period)-3)

        print("- Spin period= {} ms".format(self.ps*1000))
        print("- Orbital period= {} days".format(self.pb))
        print("- Projected axis= {} ls".format(self.x))
        print("- CHI2R= {}".format(self.chi2r))

    def fit_ellipse_circular(self):

        print("")
        print("Fitting a circular orbit (OM at 0, and T0=TA).")

        par = least_squares(residuals_circular, x0=[self.ps,self.pb,self.x], args=(self.period,self.derivative,self.der_error))

        (self.ps,self.pb,self.x)=(par.x[0],par.x[1],par.x[2])
#        self.acc_amp=np.square(2*np.pi/self.pb)*self.x*299792458  #Use the ones from the coarse fit
#        self.ps_amp=2*np.pi*self.ps*self.x/self.pb                #I mean we also do that for the eccentric fit
        residuals=residuals_circular([self.ps,self.pb,self.x],self.period,self.derivative,self.der_error)

        chi2=np.sum(residuals**2)
        self.chi2r=chi2/(np.size(self.period)-3)

        print("Estimating uncertainties.")
        print("")

        err=[0,0,0]
        par2=[self.ps,self.pb,self.x]

        for i, parameter in enumerate(par.x):
            if i==0:
                var=np.arange(parameter-self.ps_amp/10,parameter+self.ps_amp/10,self.ps_amp/10000)
            if i==1:
                dpb=self.acc_amp/(np.square(2*np.pi/86400)*self.pb**(-3/2)*self.x*299792458/2)
                var=np.arange(parameter-dpb/10,parameter+dpb/10,dpb/10000)
            if i==2:
                dx=self.acc_amp/(np.square(2*np.pi/(self.pb*86400))*299792458)
                var=np.arange(parameter-dx/10,parameter+dx/10,dx/10000)
            chi2_var=[]
            for element in var:
                par2[i]=element
                residuals=residuals_circular(par2,self.period,self.derivative,self.der_error)
                chi2_var.append(np.sum(residuals**2))
            chi2_var=np.array(chi2_var)
            mask = ~np.isnan(chi2_var)
            var = var[mask]
            chi2_var = chi2_var[mask]
            var=var[chi2_var<(np.min(chi2_var)+3)]
            chi2_var=chi2_var[chi2_var<(np.min(chi2_var)+3)]
            err[i]=abs(var[0]-var[-1])/2
            par2=[self.ps,self.pb,self.x]

        print("- Spin period= {}+/-{} ms".format(self.ps*1000,err[0]*1000))
        print("- Orbital period= {}+/-{} days".format(self.pb,err[1]))
        print("- Projected axis= {}+/-{} ls".format(self.x,err[2]))
        print("- CHI2R= {}".format(self.chi2r))

    def fit_ellipse_eccentric(self):

        print("")
        print("Fitting an eccentric orbit.")
        print("")

        ps_amp=self.ps_amp
        dpb=self.acc_amp/(np.square(2*np.pi/86400)*self.pb**(-3/2)*self.x*299792458/2)
        dx=self.acc_amp/(np.square(2*np.pi/(self.pb*86400))*299792458)

        print("Searching in the range of:")
        print("- Spin period= {} <--> {} ms".format(1000*(self.ps-3*ps_amp),1000*(self.ps+3*ps_amp)))
        print("- Orbital period= {} <--> {} days".format(self.pb-3*dpb,self.pb+6*dpb))
        print("- Projected axis= {} <--> {} ls".format(self.x-dx,self.x+3*dx))
        print("- Eccentricity= {} <--> {}".format(0.0,1.0))
        print("- Longitude of periastron= {} <--> {} deg".format(0,360))

        bounds=[(self.ps-2*ps_amp,self.ps+2*ps_amp),(self.pb-dpb,self.pb+3*dpb),(self.x-dx,self.x+3*dx),(0.0,1.0),(0,360)]

        par = differential_evolution(residuals_eccentric, bounds, args=(self.period,self.derivative,self.der_error))

        (self.ps,self.pb,self.x,self.ecc,self.om)=(par.x[0],par.x[1],abs(par.x[2]),abs(par.x[3]),par.x[4]%360)
#        self.acc_amp=np.square(2*np.pi/self.pb)*self.x*299792458    #Yes, maintain the ones from the circular fit, they may be more stable
#        self.ps_amp=2*np.pi*self.ps*self.x/self.pb                  #Either way, the formula is not correct for eccentric orbits
        chi2=residuals_eccentric([self.ps,self.pb,self.x,self.ecc,self.om],self.period,self.derivative,self.der_error)
        self.chi2r=chi2/(np.size(self.period)-5)

        print("")
        print("Estimating uncertainties.")
        print("")

        err=[0,0,0,0,0]
        par2=[self.ps,self.pb,self.x,self.ecc,self.om]

        for i, parameter in enumerate(par.x):
            if i==0:
                var=np.arange(parameter-ps_amp/5,parameter+ps_amp/5,ps_amp/10000)
            if i==1:
                dpb=self.acc_amp/(np.square(2*np.pi/86400)*self.pb**(-3/2)*self.x*299792458/2)
                var=np.arange(parameter-dpb/5,parameter+dpb/5,dpb/10000)
            if i==2:
                dx=self.acc_amp/(np.square(2*np.pi/(self.pb*86400))*299792458)
                var=np.arange(parameter-dx/5,parameter+dx/5,dx/10000)
            if i==3:
                var=np.arange(max(0,self.ecc-0.1),min(1,self.ecc+0.1),1/10000)
            if i==4:
                var=np.arange(self.om-90,self.om+90,0.01)
            chi2_var=[]
            for element in var:
                par2[i]=element
                residuals=residuals_eccentric(par2,self.period,self.derivative,self.der_error)
                chi2_var.append(np.sum(residuals**2))
            chi2_var=np.array(chi2_var)
            mask = ~np.isnan(chi2_var)
            var = var[mask]
            chi2_var = chi2_var[mask]
            var=var[chi2_var<(np.min(chi2_var)+3)]
            chi2_var=chi2_var[chi2_var<(np.min(chi2_var)+3)]
            err[i]=abs(var[0]-var[-1])/2
            par2=[self.ps,self.pb,self.x,self.ecc,self.om]

        if round(self.ecc,4)==0:
            err[3]=2*err[3]

        print("- Spin period= {}+/-{} ms".format(self.ps*1000,err[0]*1000))
        print("- Orbital period= {}+/-{} days".format(self.pb,err[1]))
        print("- Projected axis= {}+/-{} ls".format(self.x,err[2]))
        print("- Eccentricity= {}+/-{}".format(self.ecc,err[3]))
        print("- Longitude of periastron= {}+/-{} deg".format(self.om,err[4]))
        print("- CHI2R= {}".format(self.chi2r))


    def plot_ellipse_circular(self,psr_name,show=True):

        ps=self.ps
        ps_amp=self.ps_amp

        ps_values=np.arange(ps-ps_amp,ps-ps_amp*0.99,ps_amp*0.01/500)
        ps_values=np.concatenate((ps_values,np.arange(ps-ps_amp*0.99,ps-ps_amp*0.95,ps_amp*0.04/500)))
        ps_values=np.concatenate((ps_values,np.arange(ps-ps_amp*0.95,ps-ps_amp*0.8,ps_amp*0.15/500)))
        ps_values=np.concatenate((ps_values,np.arange(ps-ps_amp*0.8,ps+ps_amp*0.8,ps_amp*1.6/500)))
        ps_values=np.concatenate((ps_values,np.arange(ps+ps_amp*0.8,ps+ps_amp*0.95,ps_amp*0.15/500)))
        ps_values=np.concatenate((ps_values,np.arange(ps+ps_amp*0.95,ps+ps_amp*0.99,ps_amp*0.04/500)))
        ps_values=np.concatenate((ps_values,np.arange(ps+ps_amp*0.99,ps+ps_amp,ps_amp*0.01/500)))

        (model_pos,model_neg)=fdot_from_f_circular(1/ps_values,self.ps,self.pb,self.x)
        
        #Transform into acceleration values.
        model_pos=-299792458*model_pos*ps_values
        model_neg=-299792458*model_neg*ps_values

        plt.ylabel("Binary LOS accel., $a_\\mathrm{LOS,b}$ (m s$^{-2}$)",size=14)
        plt.xlabel(r"$P_\mathrm{bary}$ (ms)",size=14)
        plt.errorbar(self.period*1000,self.acceleration,yerr=self.acc_error,fmt="bo",label="data",markersize=5)
        plt.plot(ps_values*1000,model_pos,"c-",label="$P_\mathrm{{s}}$ = {} ms, $x$ = {} ls,\n$P_\mathrm{{b}}$ = {} days".format(round(self.ps*1000,4),round(self.x,4),round(self.pb,4)))
        plt.plot(ps_values*1000,model_neg,"c-")
        plt.legend(fontsize=12)
        plt.title(psr_name,size=14)
        plt.tick_params(labelsize=12)
        plt.tight_layout()
        if show:
            plt.show()

    def plot_ellipse_eccentric(self,psr_name,show=True):

        ps_values=np.arange(self.ps-2*self.ps_amp,self.ps+2*self.ps_amp,self.ps_amp/10000)

        #Plot.
        (model1,model2)=fdot_from_f_eccentric(1/ps_values,self.ps,self.pb,self.x,self.ecc,self.om)

        #Transform into acceleration values.
        model1=-299792458*model1*ps_values
        model2=-299792458*model2*ps_values

        plt.ylabel("Binary LOS accel., $a_\\mathrm{LOS,b}$ (m s$^{-2}$)",size=14)
        plt.xlabel(r"$P_\mathrm{SSB}$ (ms)",size=14)
        plt.errorbar(self.period*1000,self.acceleration,yerr=self.acc_error,fmt="bo",label="data",markersize=5)
        plt.plot(ps_values*1000,model1,"c-",label="$P_\mathrm{{s}}$ = {} ms, $x$ = {} ls,\n$P_\mathrm{{b}}$ = {} days, $e$ = {},\n$\omega$ = {} deg".format(round(self.ps*1000,4),round(self.x,4),round(self.pb,4),round(self.ecc,4),round(self.om,4)))
        plt.plot(ps_values*1000,model2,"c-")
        plt.legend(fontsize=12)
        plt.title(psr_name,size=14)
        plt.tick_params(labelsize=12)
        plt.tight_layout()
        if show:
            plt.show()