import numpy as np
import pymc
import sys
import pyfits
from pylens import pylens, SBModels, MassModels, plotting_tools, convolve
import pickle
from scipy.stats import truncnorm
from scipy.interpolate import splrep, splev
import emcee
from scipy.optimize import nnls


configfile = sys.argv[1]
config = pylens.read_config(configfile)

images = {}
sigmas = {}
psfs = {}
light_models = {}
source_models = {}
lens_models = []
zp = {}
pars = []
bounds = []
steps = []
filters = [filt for filt in config['filters']]
fitbands = [band for band in config['fitbands']]
rgbbands = [band for band in config['rgbbands']]

#defines model parameters
par2index = {}
index2par = {}

nlight = len(config['light_components'])
nsource = len(config['source_components'])
nlens = len(config['lens_components'])

ncomp = 0
npar = 0
for comp in config['light_components']:
    ncomp += 1
    for par in comp['pars']:
        parpar = comp['pars'][par]
        if parpar['link'] is None and parpar['var'] == 1:
            pars.append(pymc.Uniform(par+str(ncomp), lower=parpar['low'], upper=parpar['up'], value=parpar['value']))
            bounds.append((parpar['low'], parpar['up']))
            steps.append(parpar['step'])
            par2index['light'+str(ncomp)+'.'+par] = npar
            index2par[npar] = 'light'+str(ncomp)+'.'+par
            npar += 1

ncomp = 0
for comp in config['source_components']:
    ncomp += 1
    for par in comp['pars']:
        parpar = comp['pars'][par]
        if parpar['link'] is None and parpar['var'] == 1:
            pars.append(pymc.Uniform(par+str(ncomp), lower=parpar['low'], upper=parpar['up'], value=parpar['value']))
            bounds.append((parpar['low'], parpar['up']))
            steps.append(parpar['step'])
            par2index['source'+str(ncomp)+'.'+par] = npar
            index2par[npar] = 'source'+str(ncomp)+'.'+par
            npar += 1

ncomp = 0
for comp in config['lens_components']:
    ncomp += 1
    for par in comp['pars']:
        parpar = comp['pars'][par]
        if parpar['link'] is None and parpar['var'] == 1:
            pars.append(pymc.Uniform(par+str(ncomp), lower=parpar['low'], upper=parpar['up'], value=parpar['value']))
            bounds.append((parpar['low'], parpar['up']))
            steps.append(parpar['step'])
            par2index['lens'+str(ncomp)+'.'+par] = npar
            index2par[npar] = 'lens'+str(ncomp)+'.'+par
            npar += 1

npars = len(pars)

i = 0

for band in config['filters']:

    zp[band] = config['zeropoints'][i]

    hdu = pyfits.open(config['data_dir']+'/'+config['filename']+'_%s'%band+config['science_tag'])[0]

    img = hdu.data.copy()
    subimg = img.copy()
    images[band] = subimg
    suberr = pyfits.open(config['data_dir']+'/'+config['filename']+'_%s'%band+config['err_tag'])[0].data.copy()

    if config['err_type'] == 'VAR':
        sigmas[band] = suberr**0.5
    elif config['err_type'] == 'SIGMA':
        sigmas[band] = suberr.copy()
    else:
        df

    psf_file = pyfits.open(config['data_dir']+'/'+config['filename']+'_%s'%band+config['psf_tag'])
    nhdu = len(psf_file)
    found = False
    n = 0
    while not found and n < nhdu:
        if psf_file[n].data is not None:
            psf = psf_file[n].data.copy()
            found = True
        else:
            n += 1
    if not found:
        df

    m = (psf[:2].mean()+psf[-2:].mean()+psf[:,:2].mean()+psf[:,-2:].mean())/4.
    psf -= m
    psf /= psf.sum()

    psfs[band] = psf

ncomp = 1
for comp in config['lens_components']:

    name = 'lens%d'%ncomp

    pars_here = {}
    for par in comp['pars']:
        if comp['pars'][par]['link'] is None:
            if comp['pars'][par]['var'] == 1:
                pars_here[par] = pars[par2index[name+'.'+par]]
            else:
                pars_here[par] = comp['pars'][par]['value']
        else:
            pars_here[par] = pars[par2index[comp['pars'][par]['link']]]

    ncomp += 1

    lens = MassModels.PowerLaw(name, pars_here)
    lens_models.append(lens)

ncomp = 1
light_pardicts = []
for comp in config['light_components']:

    name = 'light%d'%ncomp

    pars_here = {}
    for par in comp['pars']:
        if comp['pars'][par]['link'] is None:
            if comp['pars'][par]['var'] == 1:
                pars_here[par] = pars[par2index[name+'.'+par]]
            else:
                pars_here[par] = comp['pars'][par]['value']
        else:
            pars_here[par] = pars[par2index[comp['pars'][par]['link']]]

    ncomp += 1

    light_pardicts.append(pars_here)

ncomp = 1
source_pardicts = []
for comp in config['source_components']:

    name = 'source%d'%ncomp

    pars_here = {}
    for par in comp['pars']:
        if comp['pars'][par]['link'] is None:
            if comp['pars'][par]['var'] == 1:
                pars_here[par] = pars[par2index[name+'.'+par]]
            else:
                pars_here[par] = comp['pars'][par]['value']
        else:
            pars_here[par] = pars[par2index[comp['pars'][par]['link']]]

    ncomp += 1

    source_pardicts.append(pars_here)

for band in config['filters']:

    light_models[band] = []
    source_models[band] = []

    for pardict in light_pardicts:

        light = SBModels.Sersic('light', pardict)
        light.convolve = convolve.convolve(images[band], psfs[band])[1]
        light_models[band].append(light)

    for pardict in source_pardicts:

        source = SBModels.Sersic('source', pardict)
        source.convolve = convolve.convolve(images[band], psfs[band])[1]
        source_models[band].append(source)

ny, nx = images[filters[0]].shape
X, Y = np.meshgrid(np.arange(1.*nx), np.arange(1.*ny))
R = ((X - nx/2)**2 + (Y - ny/2)**2)**0.5

if config['maskname'] is not None:
    MASK = pyfits.open(config['data_dir']+config['maskname'])[0].data.copy()
else:
    MASK = np.ones(X.shape, dtype=int)

if config['rmax'] is not None:
    MASK[R > config['rmax']] = 0

mask = MASK > 0
mask_r = mask.ravel()

output = {}

if config['do_fit'] == 'YES':
    start = []
    for j in range(npars):
        a, b = (bounds[j][0] - pars[j].value)/steps[j], (bounds[j][1] - pars[j].value)/steps[j]
        tmp = truncnorm.rvs(a, b, size=config['Nwalkers'])*steps[j] + pars[j].value

        start.append(tmp)

    start = np.array(start).T

    npars = len(pars)

    def logprior(allpars):
        for i in range(npars):
            if allpars[i] < bounds[i][0] or allpars[i] > bounds[i][1]:
                return -np.inf
        return 0.

    nwalkers = len(start)

    fakemags = {}
    for band in fitbands:
        maglist = []
        for light in light_models[band]:
            maglist.append(99.)
        for source in source_models[band]:
            maglist.append(99.)
        fakemags[band] = maglist

    def logpfunc(allpars):
        lp = logprior(allpars)
        if not np.isfinite(lp):
            return -np.inf, fakemags

        for j in range(0, npars):
            pars[j].value = allpars[j]
        sumlogp = 0.

        for lens in lens_models:
            lens.setPars()

        xl, yl = pylens.getDeflections(lens_models, (X, Y))

        mags = {}
        for band in fitbands:

            modlist = []
            maglist = []

            for light in light_models[band]:
                light.setPars()
                light.amp = 1.
                lmodel = convolve.convolve(light.pixeval(X, Y), light.convolve, False)[0]
                modlist.append((lmodel/sigmas[band]).ravel()[mask_r])

            for source in source_models[band]:
                source.setPars()
                source.amp = 1.
                smodel = convolve.convolve(source.pixeval(xl, yl), source.convolve, False)[0]
                modlist.append((smodel/sigmas[band]).ravel()[mask_r])

            modarr = np.array(modlist).T

            if np.isnan(modarr).any():

                return -1e300, fakemags

            amps, chi = nnls(modarr, (images[band]/sigmas[band]).ravel()[mask_r])

            i = 0
            for light in light_models[band]:
                if amps[i] > 0.:
                    light.amp *= amps[i]
                    maglist.append(light.Mag(zp[band]))
                else:
                    maglist.append(99.)
                i += 1
            for source in source_models[band]:
                if amps[i] > 0.:
                    source.amp *= amps[i]
                    maglist.append(source.Mag(zp[band]))
                else:
                    maglist.append(99.)
                i += 1

            mags[band] = maglist

            logp = -0.5*chi

            if logp != logp:
                return -np.inf, fakemags
            sumlogp += logp

        return sumlogp, mags

    sampler = emcee.EnsembleSampler(nwalkers, npars, logpfunc)

    print "fitting model..."

    sampler.run_mcmc(start, config['Nsteps'])

    chain = sampler.chain
    magschain = sampler.blobs

    ML = sampler.flatlnprobability.argmax()

    for j in range(0, npars):
        pars[j].value = sampler.flatchain[ML, j]

    outchain = {}
    for i in range(npars):
        outchain[index2par[i]] = chain[:, :, i]

    for i in range(nlight):
        for band in fitbands:
            outchain['light%d.mag_%s'%(i+1, band)] = np.zeros((nwalkers, config['Nsteps']))
    for i in range(nsource):
        for band in fitbands:
            outchain['source%d.mag_%s'%(i+1, band)] = np.zeros((nwalkers, config['Nsteps']))

    for i in range(config['Nsteps']):
        for j in range(nwalkers):
            for band in fitbands:
                for l in range(nlight):
                    outchain['light%d.mag_%s'%(l+1, band)][j, i] = magschain[i][j][band][l]
                for s in range(nsource):
                    outchain['source%d.mag_%s'%(s+1, band)][j, i] = magschain[i][j][band][nlight+s]

    output['chain'] = outchain

light_ml_model = {}
source_ml_model = {}

# saves best fit model images
for lens in lens_models:
    lens.setPars()

xl, yl = pylens.getDeflections(lens_models, (X, Y))

sumlogp = 0.
mags = {}
for band in filters:

    light_ml_model[band] = []
    source_ml_model[band] = []

    modlist = []

    mags[band] = []

    for light in light_models[band]:
        light.setPars()
        light.amp = 1.
        lmodel = convolve.convolve(light.pixeval(X, Y), light.convolve, False)[0]
        modlist.append((lmodel/sigmas[band]).ravel()[mask_r])

    for source in source_models[band]:
        source.setPars()
        source.amp = 1.
        smodel = convolve.convolve(source.pixeval(xl, yl), source.convolve, False)[0]
        modlist.append((smodel/sigmas[band]).ravel()[mask_r])

    modarr = np.array(modlist).T

    amps, chi = nnls(modarr, (images[band]/sigmas[band]).ravel()[mask_r])

    if band in fitbands:
        sumlogp += -0.5*chi

    n = 0
    for light in light_models[band]:
        lmodel = convolve.convolve(light.pixeval(X, Y), light.convolve, False)[0]
        lmodel *= amps[n]
        light.amp *= amps[n]
        
        mags[band].append(light.Mag(zp[band]))
        
        light_ml_model[band].append(lmodel)
        n += 1

    for source in source_models[band]:
        smodel = convolve.convolve(source.pixeval(xl, yl), source.convolve, False)[0]
        smodel *= amps[n]
        if amps[n] > 0.:
            source.amp *= amps[n]
            mags[band].append(source.Mag(zp[band]))
        else:
            mags[band].append(99.)

        source_ml_model[band].append(smodel)
        n += 1

# makes model rgb
if len(rgbbands) >= 3:
    bandshere = rgbbands
else:
    bandshere = filters
    if len(bandshere) < 3:
        bandshere += filters
    if len(bandshere) < 3:
        bandshere += filters

sci_list = []
light_list = []
source_list = []
for band in bandshere:
    sci_list.append(images[band])
    lmodel = 0.*images[band]
    smodel = 0.*images[band]
    for light in light_ml_model[band]:
        lmodel += light
    light_list.append(lmodel)

    for source in source_ml_model[band]:
        smodel += source
    source_list.append(smodel)

plotting_tools.make_model_rgb(sci_list, light_list, source_list, outname=config['output_dir']+'/'+config['rgbname'])

output['light_ml_model'] = light_ml_model
output['source_ml_model'] = source_ml_model
output['logp'] = sumlogp

f = open(config['output_dir']+'/'+config['outname'], 'w')
pickle.dump(output, f)
f.close()

# writes a new configuration file
conflines = []
confpars = ['data_dir', 'output_dir', 'filename', 'science_tag', 'err_tag', 'err_type', 'psf_tag', 'rmax', 'Nwalkers', 'Nsteps']
for parname in confpars:
    if config[parname] is not None:
        conflines.append('%s: %s\n'%(parname, config[parname]))
filtline = 'filters: '
zpline = 'zeropoints: '
nfilt = len(filters)
for i in range(nfilt-1):
    filtline += '%s, '%filters[i]
    zpline += '%f, '%zp[filters[i]]
filtline += filters[-1]+'\n'
zpline += '%f\n'%zp[filters[-1]]
conflines.append(filtline)
conflines.append(zpline)

filtline = 'fitbands: '
nfilt = len(fitbands)
for i in range(nfilt-1):
    filtline += '%s, '%fitbands[i]
filtline += fitbands[-1]+'\n'
conflines.append(filtline)

if config['rgbbands'] is not None:
    filtline = 'rgbbands: '
    nfilt = len(rgbbands)
    for i in range(nfilt-1):
        filtline += '%s, '%rgbbands[i]
    filtline += rgbbands[-1]+'\n'
    conflines.append(filtline)

conflines.append('\n')
conflines.append('# MODELS\n')

serpars = ['x', 'y', 'pa', 'q', 're', 'n']
ncomp = 0
for light in light_pardicts:
    conflines.append('\n')
    conflines.append('light_model Sersic\n')
    for par in serpars:
        parname = 'light%d.%s'%(ncomp+1, par)
        if parname in par2index:
            npar = par2index[parname]
            conflines.append('%s %f %f %f %f 1\n'%(par, pars[npar].value, bounds[npar][0], bounds[npar][1], steps[npar]))
        else:
            if config['light_components'][ncomp]['pars'][par]['link'] is None:
                conflines.append('%s %f -1 -1 -1 0\n'%(par, config['light_components'][ncomp]['pars'][par]['value']))
            else:
                lname = config['light_components'][ncomp]['pars'][par]['link']
                npar = par2index[lname]
                conflines.append('%s %f %f %f %f 1 %s\n'%(par, pars[npar].value, bounds[npar][0], bounds[npar][1], steps[npar], lname))
    for band in filters:
        conflines.append('mag_%s %3.2f\n'%(band, mags[band][ncomp]))
    ncomp += 1

ncomp = 0
for source in source_pardicts:
    conflines.append('\n')
    conflines.append('source_model Sersic\n')
    for par in serpars:
        parname = 'source%d.%s'%(ncomp+1, par)
        if parname in par2index:
            npar = par2index[parname]
            conflines.append('%s %f %f %f %f 1\n'%(par, pars[npar].value, bounds[npar][0], bounds[npar][1], steps[npar]))
        else:
            if config['source_components'][ncomp]['pars'][par]['link'] is None:
                conflines.append('%s %f -1 -1 -1 0\n'%(par, config['source_components'][ncomp]['pars'][par]['value']))
            else:
                lname = config['source_components'][ncomp]['pars'][par]['link']
                npar = par2index[lname]
                conflines.append('%s %f %f %f %f 1 %s\n'%(par, pars[npar].value, bounds[npar][0], bounds[npar][1], steps[npar], lname))
    for band in filters:
        conflines.append('mag_%s %3.2f\n'%(band, mags[band][nlight+ncomp]))
    ncomp += 1

ncomp = 0
powpars = ['x', 'y', 'pa', 'q', 'b', 'eta']

for lens in source_pardicts:
    conflines.append('\n')
    conflines.append('lens_model Powerlaw\n')
    for par in powpars:
        parname = 'lens%d.%s'%(ncomp+1, par)
        if parname in par2index:
            npar = par2index[parname]
            conflines.append('%s %f %f %f %f 1\n'%(par, pars[npar].value, bounds[npar][0], bounds[npar][1], steps[npar]))
        else:
            if config['lens_components'][ncomp]['pars'][par]['link'] is None:
                conflines.append('%s %f -1 -1 -1 0\n'%(par, config['lens_components'][ncomp]['pars'][par]['value']))
            else:
                lname = config['lens_components'][ncomp]['pars'][par]['link']
                npar = par2index[lname]
                conflines.append('%s %f %f %f %f 1 %s\n'%(par, pars[npar].value, bounds[npar][0], bounds[npar][1], steps[npar], lname))
    ncomp += 1

conflines.append('\n')
conflines.append('logp %f\n'%sumlogp)

f = open(config['output_dir']+'/'+configfile+'.out', 'w')
f.writelines(conflines)
f.close()

