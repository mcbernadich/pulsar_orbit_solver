# pulsar_orbit_solver

Some tools to solve puslar orbits before timing. These work on pdmp or perpfold measurements of spin periods, spin period derivatives, and accelerataions.

## Included scripts

```orbital_tools.py```: a library with classes and functions. See headers for a descritpion.
```fitOrbit.py```: an interactive Keplerian orbit fitter that reads in parameter files and MJD - spin data. Usually requires a good estimate of the orbital period to work.
```searchOrbit.py```: a script that searches for the best orbital period over a selected range from MJD - spin data. It implements Keplerian fits on a grid of Pb values. Good to find the orbital period to use as a prior for ```fitOrbit.py```
```ellipseFit.py```: a script that fits Keplerian aprameters to a set of spin - acceleration measurements ([Freire et al. 2001](https://ui.adsabs.harvard.edu/abs/2001MNRAS.322..885F/abstract)). Good to estimate parameters to use as priors in ```searchOrbit.py``` and/or ```fitOrbit.py```.

Try ```python3 ${script} -h``` to see the option

## Example data

Some example data is included to test these scripts

### J1326-4728K from Omega Centauri

This includes a prior parameter file ```J1326-4728K.par``` and, MJD, Ps (ms), dPs (ms) measurements from individual observations in ```*.per``` files (from pdmp). It has a a circular orbit, so you can test fitting circular orbits with. The quality of the measurements can be very low, however, so expect some large reduced chi^2 values.

Try ```python3 searchOrbit.py -d "J1326-4728K/*.per" -p J1326-4728K/J1326-4728K.par -r "0.08:0.10"``` to find its orbital period.

Try ```python3 fitOrbit.py -d "J1326-4728K/*.per" -p J1326-4728K/J1326-4728K.par``` to make a full Keplerian fit.

### J1326-4728N from Omega Centauri

This includes a prior parameter file ```J1326-4728N.par``` and, MJD, Ps (ms), dPs (ms) measurements from individual observations in ```*.per``` files (from pdmp). In addition, it also includes a ```J1326-4728K.dat``` file with MJD, F0, F1 measurements (from tempo2 on individual osbervations). This pulsar has an eccentric orbit, so you will need to include eccentricity in your fits.

Try ```python3 fitEllipse.py --freq 1 -d J1326-4728N/J1326-4728N.dat``` to get some preliminary orbital parameters and visualise its eccentricity.

Try ```python3 searchOrbit.py -d "J1326-4728N/*.per" -p J1326-4728N/J1326-4728N.par -r "5:8" -e "0.05:0.15"``` to search for its best orbital period. Notice that the Pb and eccentricity ranges are derived from the results of the previous fit. The example uses ```J1326-4728N.par```, but you could make your own parameter file from the result of the previous fit.

Try ```python3 fitOrbit.py -d "J1326-4728N/*.per" -p J1326-4728N/J1326-4728N.par``` to perform the full Keplerian fit. Notice that, again, ```J1326-4728N.par``` is used as a prior, but you could alos use your own prior from the previous two fits.

## Feedback

Write to mcbernadich@mpifr-bonn.mpg.de or miquel.colomibernadich@inaf.it for feeddback

## Credits

Soon a paper will be submited where this code is used