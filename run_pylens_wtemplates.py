import numpy as np
import sys
import pyfits
from pylens import pylens, SBModels, MassModels, plotting_tools, convolve
import pickle
from scipy.stats import truncnorm
from scipy.interpolate import splrep, splev
import emcee
from scipy.optimize import nnls
import os


configfile = sys.argv[1]
rootdir = os.environ.get('PYLENSDIR')

class Par:

    def __init__(self, name, lower=0., upper=1., value=0.):

        self.name = name
        self.lower = lower
        self.upper = upper
        self.value = value

def read_config(filename):

    f = open(filename, 'r')
    lines = f.readlines()
    f.close()

    config = {'data_dir':'./', 'mask_dir': None, 'output_dir':'./', 'filters': None, 'main_band': None, 'fitbands': None, 'rgbbands': None, \
              'zeropoints': None, \
              'filename': None, 'filter_prefix': '', 'filter_suffix': '', 'science_tag':'_sci.fits', 'err_tag':'_var.fits', 'err_type': 'VAR', 'psf_tag':'_psf.fits', \
              'rmax': None, 'do_fit': 'YES', 'Nsteps':300, 'Nwalkers':30, 'burnin':None, 'maskname':None, \
              'rgbname': None, 'rgbcuts': None, 'outname': None}

    preamble = True

    i = 0
    while preamble and i < len(lines):
        if '#' in lines[i] and 'MODELS' in lines[i]:
            preamble = False
        else:
            line = lines[i].split('#')[0].split()
            if len(line) > 0:
                parname = line[0].split(':')[0]
                if parname in config:
                    config[parname] = lines[i].split('#')[0].split(':')[1].split('\n')[0].lstrip().rstrip()
        i += 1

    filtlist = []
    filternames = config['filters'].split(',')
    for name in filternames:
        filtlist.append(name.lstrip())
    config['filters'] = filtlist
    if config['fitbands'] is not None:
        filtlist = []
        filternames = config['fitbands'].split(',')
        for name in filternames:
            filtlist.append(name.lstrip())
        config['fitbands'] = filtlist
    else:
        config['fitbands'] = config['filters']

    config['colors'] = []
    for band in config['fitbands']:
        if band != config['main_band']:
            config['colors'].append('%s-%s'%(band, config['main_band']))

    if config['rgbcuts'] is not None:
        cutlist = []
        cuts = config['rgbcuts'].split(',')
        for cut in cuts:
            cutlist.append(float(cut))

        config['rgbcuts'] = cutlist
    else:
        config['rgbcuts'] = (99., 99., 99.)

    if config['outname'] is None:
        config['outname'] = filename+'.output'

    rgblist = []
    rgbnames = config['rgbbands'].split(',')
    for name in rgbnames:
        rgblist.append(name.lstrip())
    config['rgbbands'] = rgblist
    config['zeropoints'] = np.array(config['zeropoints'].split(','), dtype='float')
    config['Nsteps'] = int(config['Nsteps'])
    config['Nwalkers'] = int(config['Nwalkers'])

    light_components = []
    lens_components = []
    source_templates = []

    while i < len(lines):

        line = lines[i].split()
        if len(line) > 0:
            if line[0] == 'light_model':
                model_class = line[1]

                if model_class == 'Sersic':
                    npars = 6 + len(config['colors'])
                    parnames = ['x', 'y', 'q', 'pa', 're', 'n']
                    parnames += config['colors']
                else:
                    df

                comp = {'class':model_class, 'pars':{}}

                foundpars = 0
                j = 1

                while foundpars < npars and j+i < len(lines):
                    line = lines[j+i].split()
                    if lines[j+i][0] != '#' and len(line) > 0:
                        if line[0] in parnames:
                            foundpars += 1
                            par = line[0]
                            link = None
                            if len(line) > 6:
                                link = line[6]
                            tmp_par = {'value': float(line[1]), 'low': float(line[2]), 'up': float(line[3]), \
                           'step': float(line[4]), 'var': int(line[5]), 'link':link}
                            comp['pars'][par] = tmp_par
                    j += 1

                i += j

                if foundpars < npars:
                    print 'not all parameters found!'
                else:
                    light_components.append(comp)

            elif 'source_template' in line[0]:
                model_class = line[1].lstrip()
                tempname = line[2].rstrip()

                if model_class == 'Sersic':
                    npars = 7
                    parnames = ['x', 'y', 'q', 'pa', 're', 'n', 'zs']
                else:
                    df

                comp = {'class':model_class, 'pars':{}, 'tempname': tempname}

                foundpars = 0
                j = 1

                while foundpars < npars and j+i < len(lines):
                    line = lines[j+i].split()
                    if lines[j+i][0] != '#' and len(line) > 0:
                        if line[0] in parnames:
                            foundpars += 1
                            par = line[0]
                            link = None
                            if len(line) > 6:
                                link = line[6]
                            tmp_par = {'value': float(line[1]), 'low': float(line[2]), 'up': float(line[3]), \
                           'step': float(line[4]), 'var': int(line[5]), 'link':link}
                            comp['pars'][par] = tmp_par
                    j += 1

                i += j

                if foundpars < npars:
                    print 'not all parameters found!'
                else:
                    source_templates.append(comp)

            elif 'lens_model' in line[0]:
                model_class = line[1].lstrip()

                if model_class == 'Powerlaw':
                    npars = 6
                    parnames = ['x', 'y', 'q', 'pa', 'b', 'eta']
                else:
                    df

                comp = {'class':model_class, 'pars':{}}

                foundpars = 0
                j = 1

                while foundpars < npars and j+i < len(lines):
                    line = lines[j+i].split()
                    if lines[j+i][0] != '#' and len(line) > 0:
                        if line[0] in parnames:
                            foundpars += 1
                            par = line[0]
                            link = None
                            if len(line) > 6:
                                link = line[6]
                            tmp_par = {'value': float(line[1]), 'low': float(line[2]), 'up': float(line[3]), \
                           'step': float(line[4]), 'var': int(line[5]), 'link':link}
                            comp['pars'][par] = tmp_par
                    j += 1

                i += j

                if foundpars < npars:
                    print 'not all parameters found!'
                else:
                    lens_components.append(comp)

            else:
                i += 1
        else:
            i += 1

    config['light_components'] = light_components
    config['source_templates'] = source_templates
    config['lens_components'] = lens_components

    return config

config = read_config(configfile)

images = {}
sigmas = {}
convol_matrix = {}
light_models = []
sourcetemp_models = []
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
nsource = len(config['source_templates'])
nlens = len(config['lens_components'])

ncomp = 0
npar = 0
for comp in config['light_components']:
    ncomp += 1
    for par in comp['pars']:
        parpar = comp['pars'][par]
        if parpar['link'] is None and parpar['var'] == 1:
            pars.append(Par(par+str(ncomp), lower=parpar['low'], upper=parpar['up'], value=parpar['value']))
            bounds.append((parpar['low'], parpar['up']))
            steps.append(parpar['step'])
            par2index['light'+str(ncomp)+'.'+par] = npar
            index2par[npar] = 'light'+str(ncomp)+'.'+par
            npar += 1

ncomp = 0
for comp in config['source_templates']:
    ncomp += 1
    for par in comp['pars']:
        parpar = comp['pars'][par]
        if parpar['link'] is None and parpar['var'] == 1:
            pars.append(Par(par+str(ncomp), lower=parpar['low'], upper=parpar['up'], value=parpar['value']))
            bounds.append((parpar['low'], parpar['up']))
            steps.append(parpar['step'])
            par2index['sourcetemp'+str(ncomp)+'.'+par] = npar
            index2par[npar] = 'sourcetemp'+str(ncomp)+'.'+par
            npar += 1

ncomp = 0
for comp in config['lens_components']:
    ncomp += 1
    for par in comp['pars']:
        parpar = comp['pars'][par]
        if parpar['link'] is None and parpar['var'] == 1:
            pars.append(Par(par+str(ncomp), lower=parpar['low'], upper=parpar['up'], value=parpar['value']))
            bounds.append((parpar['low'], parpar['up']))
            steps.append(parpar['step'])
            par2index['lens'+str(ncomp)+'.'+par] = npar
            index2par[npar] = 'lens'+str(ncomp)+'.'+par
            npar += 1

npars = len(pars)

i = 0

filtdic = {}
for band in config['filters']:

    zp[band] = config['zeropoints'][i]

    filtname = rootdir+'pylens/filters/%s%s%s'%(config['filter_prefix'], band, config['filter_suffix'])

    f = open(filtname, 'r')
    ftable = np.loadtxt(f)
    f.close()

    filtdic[band] = (ftable[:, 0], ftable[:, 1])

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

    convol_matrix[band] = convolve.convolve(images[band], psf)[1]

ncomp = 1
for comp in config['lens_components']:

    name = 'lens%d'%ncomp

    pars_here = {}
    for par in comp['pars']:
        if comp['pars'][par]['link'] is None:
            if comp['pars'][par]['var'] == 1:
                pars_here[par] = pars[par2index[name+'.'+par]]
            elif comp['pars'][par]['var'] == 0:
                pars_here[par] = Par(par, lower=comp['pars'][par]['value'], upper=comp['pars'][par]['value'], value=comp['pars'][par]['value'])
            else:
                df
        else:
            pars_here[par] = pars[par2index[comp['pars'][par]['link']]]

    ncomp += 1

    lens = MassModels.PowerLaw(name, pars_here)
    lens_models.append(lens)

ncomp = 1
light_pardicts = []
light_colordicts = []
sersicpars = ['x', 'y', 'q', 're', 'pa', 'n']
for comp in config['light_components']:

    name = 'light%d'%ncomp

    pars_here = {}
    colors_here = {}
    for par in comp['pars']:
        if par in sersicpars:
            if comp['pars'][par]['link'] is None:
                if comp['pars'][par]['var'] == 1:
                    pars_here[par] = pars[par2index[name+'.'+par]]
                else:
                    pars_here[par] = Par(par, lower=comp['pars'][par]['value'], upper=comp['pars'][par]['value'], value=comp['pars'][par]['value'])
            else:
                pars_here[par] = pars[par2index[comp['pars'][par]['link']]]
        else:
            if comp['pars'][par]['link'] is None:
                if comp['pars'][par]['var'] == 1:
                    colors_here[par] = pars[par2index[name+'.'+par]]
                else:
                    colors_here[par] = Par(par, lower=comp['pars'][par]['value'], upper=comp['pars'][par]['value'], value=comp['pars'][par]['value'])
            else:
                colors_here[par] = pars[par2index[comp['pars'][par]['link']]]

    ncomp += 1

    light_pardicts.append(pars_here)
    light_colordicts.append(colors_here)

ncomp = 1
sourcetemp_pardicts = []
sourcetemp_templates = []
for comp in config['source_templates']:

    name = 'sourcetemp%d'%ncomp

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

    sourcetemp_pardicts.append(pars_here)
    sourcetemp_templates.append(comp['tempname'])

for pardict, tempname in zip(sourcetemp_pardicts, sourcetemp_templates):

    f = open(rootdir+'pylens/templates/'+tempname, 'r')
    ttable = np.loadtxt(f)
    f.close()

    template = (ttable[:, 0], ttable[:, 1])

    source = SBModels.SersicTemplate('source', pardict, template, filtdic)
    sourcetemp_models.append(source)

for pardict in light_pardicts:
    light = SBModels.Sersic('light', pardict)
    light_models.append(light)

ny, nx = images[filters[0]].shape
X, Y = np.meshgrid(np.arange(1.*nx), np.arange(1.*ny))
R = ((X - nx/2)**2 + (Y - ny/2)**2)**0.5

if config['mask_dir'] is None:
    config['mask_dir'] = config['data_dir']
if config['maskname'] is not None:
    MASK = pyfits.open(config['mask_dir']+config['maskname'])[0].data.copy()
else:
    MASK = np.ones(X.shape, dtype=int)

if config['rmax'] is not None:
    MASK[R > config['rmax']] = 0

mask = MASK > 0
mask_r = mask.ravel()

output = {}

nfitbands = len(fitbands)

modelstack = np.zeros((nfitbands * ny, nx))
datastack = 0. * modelstack
sigmastack = 0. * modelstack
maskstack = np.zeros((nfitbands * ny, nx), dtype=bool)
for i in range(nfitbands):
    datastack[i*ny: (i+1)*ny, :] = images[fitbands[i]]
    sigmastack[i*ny: (i+1)*ny, :] = sigmas[fitbands[i]]
    maskstack[i*ny: (i+1)*ny, :] = mask
maskstack_r = maskstack.ravel()

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

    fakemags = []
    for comp in light_models + sourcetemp_models:
        magdic = {}
        for band in fitbands:
            magdic[band] = 99.
        fakemags.append(magdic)

    ncomp = len(light_models) + len(sourcetemp_models)

    def logpfunc(allpars):
        lp = logprior(allpars)
        if not np.isfinite(lp):
            return -np.inf, fakemags

        for j in range(0, npars):
            pars[j].value = allpars[j]

        for lens in lens_models:
            lens.setPars()

        xl, yl = pylens.getDeflections(lens_models, (X, Y))

        mags = {}

        modlist = []

        for light, colordict in zip(light_models, light_colordicts):
            lmodel = 0. * datastack

            light.setPars()
            light.amp = 1.
            lpix = light.pixeval(X, Y)
            for i in range(nfitbands):
                if fitbands[i] == config['main_band']:
                    scale = 1.
                else:
                    color = colordict['%s-%s'%(fitbands[i], config['main_band'])].value
                    scale = 10.**(-(color + zp[config['main_band']] - zp[fitbands[i]])/2.5)
                    
                lmodel[i*ny: (i+1)*ny, :] = scale * convolve.convolve(lpix, convol_matrix[fitbands[i]], False)[0]
            modlist.append((lmodel/sigmastack).ravel()[maskstack_r])

        for source in sourcetemp_models:
            smodel = 0. * datastack
            source.setPars()
            source.amp = 1.
            spix = source.pixeval(xl, yl)
            for i in range(nfitbands):
                if fitbands[i] == config['main_band']:
                    scale = 1.
                else:
                    scale = source.scale(fitbands[i], config['main_band']) * 10.**(-(zp[config['main_band']] - zp[fitbands[i]])/2.5)
                smodel[i*ny: (i+1)*ny, :] = scale * convolve.convolve(spix, convol_matrix[fitbands[i]], False)[0]
            modlist.append((smodel/sigmastack).ravel()[maskstack_r])

        modarr = np.array(modlist).T

        if np.isnan(modarr).any():

            return -1e300, fakemags

        amps, chi = nnls(modarr, (datastack/sigmastack).ravel()[maskstack_r])

        maglist = []
        i = 0
        for comp, colordict in zip(light_models, light_colordicts):
            magdic = {}
            comp.amp *= amps[i]
            for band in fitbands:
                if band == config['main_band']:
                    magdic[band] = comp.Mag(zp[config['main_band']])
                else:
                    magdic[band] = comp.Mag(zp[config['main_band']]) + colordict['%s-%s'%(band, config['main_band'])].value
            maglist.append(magdic)
            i += 1

        for comp in sourcetemp_models:
            magdic = {}
            if amps[i] > 0.:
                comp.amp *= amps[i]
                mainmag = comp.Mag(zp[config['main_band']])
                for band in fitbands:
                    if band == config['main_band']:
                        magdic[band] = mainmag
                    else:
                        scale = comp.scale(band, config['main_band']) * 10.**(-(zp[config['main_band']] - zp[fitbands[i]])/2.5)
                        magdic[band] = mainmag - 2.5*np.log10(scale)
            else:
                for band in fitbands:
                    magdic[band] = 99.
            maglist.append(magdic)
            i += 1

        logp = -0.5*chi
        if logp != logp:
            return -np.inf, fakemags

        return logp, maglist

    sampler = emcee.EnsembleSampler(nwalkers, npars, logpfunc)

    print "fitting model..."

    sampler.run_mcmc(start, config['Nsteps'])

    chain = sampler.chain
    magschain = sampler.blobs

    ML = sampler.flatlnprobability.argmax()

    for j in range(0, npars):
        pars[j].value = sampler.flatchain[ML, j]

    outchain = {}
    outchain['logp'] = sampler.lnprobability

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
                    outchain['light%d.mag_%s'%(l+1, band)][j, i] = magschain[i][j][l][band]
                for s in range(nsource):
                    outchain['source%d.mag_%s'%(s+1, band)][j, i] = magschain[i][j][nlight+s][band]

    output['chain'] = outchain

light_ml_model = []
source_ml_model = []
light_mags = []
source_mags = []

# saves best fit model images
for lens in lens_models:
    lens.setPars()

xl, yl = pylens.getDeflections(lens_models, (X, Y))

modlist = []

ntotbands = len(filters)

modelstack = np.zeros((ntotbands * ny, nx))
datastack = 0. * modelstack
sigmastack = 0. * modelstack
maskstack = np.zeros((ntotbands * ny, nx), dtype=bool)

for i in range(ntotbands):
    datastack[i*ny: (i+1)*ny, :] = images[filters[i]]
    sigmastack[i*ny: (i+1)*ny, :] = sigmas[filters[i]]
    maskstack[i*ny: (i+1)*ny, :] = mask

maskstack_r = maskstack.ravel()

for light, colordict in zip(light_models, light_colordicts):
    light.setPars()
    light.amp = 1.
    lpix = light.pixeval(X, Y)
    lmodel = 0. * datastack
    for i in range(ntotbands):
        if filters[i] == config['main_band']:
            scale = 1.
        else:
            color = colordict['%s-%s'%(filters[i], config['main_band'])].value
            scale = 10.**(-(color + zp[config['main_band']] - zp[filters[i]])/2.5)
        lmodel[i*ny: (i+1)*ny, :] = scale * convolve.convolve(lpix, convol_matrix[filters[i]], False)[0]
    modlist.append((lmodel/sigmastack).ravel()[maskstack_r])

for source in sourcetemp_models:
    source.setPars()
    source.amp = 1.
    spix = source.pixeval(xl, yl)
    smodel = 0. * datastack
    for i in range(ntotbands):
        if fitbands[i] == config['main_band']:
            scale = 1.
        else:
            scale = source.scale(filters[i], config['main_band']) * 10.**(-(zp[config['main_band']] - zp[fitbands[i]])/2.5)
        smodel[i*ny: (i+1)*ny, :] = scale * convolve.convolve(spix, convol_matrix[filters[i]], False)[0]
    modlist.append((smodel/sigmastack).ravel()[maskstack_r])

modarr = np.array(modlist).T

amps, chi = nnls(modarr, (datastack/sigmastack).ravel()[maskstack_r])

n = 0
for light, colordict in zip(light_models, light_colordicts):

    light_ml_dic = {}
    light_ml_mag = {}

    light.amp = amps[n]
    lpix = light.pixeval(X, Y)
    for band in filters:
        if band == config['main_band']:
            scale = 1.
        else:
            color = colordict['%s-%s'%(band, config['main_band'])].value
            scale = 10.**(-(color + zp[config['main_band']] - zp[band])/2.5)
        light_ml_dic[band] = scale * convolve.convolve(lpix, convol_matrix[band], False)[0]
        if band == config['main_band']:
            light_ml_mag[band] = light.Mag(zp[band])
        else:
            light_ml_mag[band] = light.Mag(zp[config['main_band']]) + colordict['%s-%s'%(band, config['main_band'])].value
        
    light_ml_model.append(light_ml_dic)
    light_mags.append(light_ml_mag)
    n += 1

for source in sourcetemp_models:

    source_ml_dic = {}
    source_ml_mag = {}

    source.amp = amps[n]
    spix = source.pixeval(xl, yl)
    mainmag = source.Mag(zp[config['main_band']])

    for band in filters:
        if band == config['main_band']:
            scale = 1.
            source_ml_mag[band] = mainmag
        else:
            scale = source.scale(band, config['main_band']) * 10.**(-(zp[config['main_band']] - zp[fitbands[i]])/2.5)
            source_ml_mag[band] = mainmag - 2.5*np.log10(scale)
        source_ml_dic[band] = scale * convolve.convolve(spix, convol_matrix[band], False)[0]

    source_ml_model.append(source_ml_dic)
    source_mags.append(source_ml_mag)
    n += 1

# makes model rgb
sci_list = []
light_list = []
source_list = []
for band in filters:
    sci_list.append(images[band])
    lmodel = 0.*images[band]
    smodel = 0.*images[band]
    for light in light_ml_model:
        lmodel += light[band]
    light_list.append(lmodel)
    for source in source_ml_model:
        smodel += source[band]
    source_list.append(smodel)

if config['rgbname'] is None:
    config['rgbname'] = config['filename'] + '_rgb.png'
plotting_tools.make_full_rgb(sci_list, light_list, source_list, outname=config['output_dir']+'/'+config['rgbname'])

"""
if len(rgbbands) >= 3:
    bandshere = rgbbands
else:
    bandshere = filters
    if len(bandshere) < 3:
        bandshere += filters
    if len(bandshere) < 3:
        bandshere += filters
rgbsets = []
if ntotbands == 1:
    rgbsets.append((filters[0], filters[0], filters[0]))
elif ntotbands == 2:
    rgbsets.append((filters[1], filters[1], filters[0]))
elif ntotbands == 3:
    rgbsets.append((filters[2], filters[1], filters[0]))
else:
    nsets = ntotbands - 2
    for i in range(nsets):
        rgbsets.append((filters[ntotbands-1-i], filters[ntotbands-2-i], filters[ntotbands-3-i]))

def make_rgb(rgbbands):
    sci_list = []
    light_list = []
    source_list = []
    for band in rgbbands:
        sci_list.append(images[band])
        lmodel = 0.*images[band]
        smodel = 0.*images[band]
        for light in light_ml_model:
            lmodel += light[band]
        light_list.append(lmodel)
        for source in source_ml_model:
            smodel += source[band]
        source_list.append(smodel)

    plotting_tools.make_model_rgb(sci_list, light_list, source_list, outname=config['output_dir']+'/'+config['outname']+'_%s%s%s.png'%rgbbands)

for rgbset in rgbsets:
    make_rgb(rgbset)

#plotting_tools.make_model_rgb(sci_list, light_list, source_list, outname=config['output_dir']+'/'+config['rgbname'])
"""

output['light_ml_model'] = light_ml_model
output['source_ml_model'] = source_ml_model
output['logp'] = -0.5*chi

f = open(config['output_dir']+'/'+config['outname'], 'w')
pickle.dump(output, f)
f.close()

# writes a new configuration file
conflines = []
confpars = ['data_dir', 'mask_dir', 'output_dir', 'filename', 'science_tag', 'err_tag', 'err_type', 'psf_tag', 'rmax', 'Nwalkers', 'Nsteps', 'main_band', 'filter_prefix', 'filter_suffix']
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

lightpars = ['x', 'y', 'pa', 'q', 're', 'n']
lightpars += config['colors']
ncomp = 0
for light, mags in zip(light_pardicts, light_mags):
    conflines.append('\n')
    conflines.append('light_model Sersic\n')
    for par in lightpars:
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
        conflines.append('mag_%s %3.2f\n'%(band, mags[band]))
    ncomp += 1

ncomp = 0
sourcepars = ['x', 'y', 'pa', 'q', 're', 'n', 'zs']
for source, mags in zip(sourcetemp_pardicts, source_mags):
    conflines.append('\n')
    conflines.append('source_model Sersic\n')
    for par in sourcepars:
        parname = 'sourcetemp%d.%s'%(ncomp+1, par)
        if parname in par2index:
            npar = par2index[parname]
            conflines.append('%s %f %f %f %f 1\n'%(par, pars[npar].value, bounds[npar][0], bounds[npar][1], steps[npar]))
        else:
            if config['source_templates'][ncomp]['pars'][par]['link'] is None:
                conflines.append('%s %f -1 -1 -1 0\n'%(par, config['source_templates'][ncomp]['pars'][par]['value']))
            else:
                lname = config['source_templates'][ncomp]['pars'][par]['link']
                npar = par2index[lname]
                conflines.append('%s %f %f %f %f 1 %s\n'%(par, pars[npar].value, bounds[npar][0], bounds[npar][1], steps[npar], lname))
    for band in filters:
        conflines.append('mag_%s %3.2f\n'%(band, mags[band]))
    ncomp += 1

ncomp = 0
powpars = ['x', 'y', 'pa', 'q', 'b', 'eta']

for lens in lens_models:
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
conflines.append('logp %f\n'%(-0.5*chi))

f = open(config['output_dir']+'/'+configfile+'.out', 'w')
f.writelines(conflines)
f.close()

