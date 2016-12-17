# MCMC-based photometry fitting and lens modeling code

Based on Matt Auger's pylens. Made a bit more user-friendly.

Can only fit Sersic profiles and singular power-law elliptical mass distribution (SPEMD) at the moment.

USAGE: python run_pylens.py configfile

If you wish you can define a $PYLENSPATH and from anywhere you can run

python $PYLENSPATH/run_pylens.py configfile


## Example

In the directory example/ there is one working example.
It fits a sersic profile to CFHT u, g, r, i, z data of a strong lens.
To run the example, from the example/ directory run

python ../run_pylens.py sersic_light

It will produce a bunch of output.

sersic_light.out is a new configuration file with the same setting as the original, except for the starting point of the model, which is set to the maximum likelihood of the MCMC chain.

sersic_light.output is a pickle'd object. To open it,

import pickle

f = open('sersic_light.output', 'r')
output = pickle.load(f)
f.close()

output is a dictionary. It contains, among other stuff, the chains for each model parameter, as well as the magnitudes in the bands used to do the fit (defined by keyword 'fitbands' in configuration file).
The MCMC is done using emcee. emcee uses a number of walkers to explore the posterior PDF (specified by keyword Nwalkers in configuration file).
In the end the chain will have shape (Nwalkers, Nsteps).


