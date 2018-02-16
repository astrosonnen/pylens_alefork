import SBProfiles as SBProfiles
from math import pi
import numpy as np
from scipy.interpolate import splrep, splev
import os


tempdir = os.environ.get('PYLENSDIR') + 'pylens/templates/'

def cnts2mag(cnts,zp):
    from math import log10
    return -2.5*log10(cnts) + zp


def etau_madau(wl, z):
    """
    Madau 1995 extinction for a galaxy spectrum at redshift z 
    defined on a wavelength grid wl
    """
    n=len(wl)
    l=np.array([1216.,1026.,973.,950.])
    xe=1.+z
    
    #If all the spectrum is redder than (1+z)*wl_lyman_alfa 
    if wl[0]> l[0]*xe: return np.zeros(n)+1.
    
    #Madau coefficients
    c=np.array([3.6e-3,1.7e-3,1.2e-3,9.3e-4])
    ll=912.
    tau=wl*0.
    i1=np.searchsorted(wl,ll)
    i2=n-1
    #Lyman series absorption
    for i in range(len(l)):
	i2=np.searchsorted(wl[i1:i2],l[i]*xe)
	tau[i1:i2]=tau[i1:i2]+c[i]*(wl[i1:i2]/l[i])**3.46

    if ll*xe < wl[0]:
        return np.exp(-tau)

    #Photoelectric absorption
    xe=1.+z
    i2=np.searchsorted(wl,ll*xe)
    xc=wl[i1:i2]/ll
    xc3=xc**3
    tau[i1:i2]=tau[i1:i2]+\
                (0.25*xc3*(xe**.46-xc**0.46)\
                 +9.4*xc**1.5*(xe**0.18-xc**0.18)\
                 -0.7*xc3*(xc**(-1.32)-xe**(-1.32))\
                 -0.023*(xe**1.68-xc**1.68))

    tau = np.clip(tau, 0, 700)
    return np.exp(-tau)
    # if tau>700. : return 0.
    # else: return exp(-tau)


class SBModel:
    def __init__(self,name,pars,convolve=0):
        if 'amp' not in pars.keys() and 'logamp' not in pars.keys():
            pars['amp'] = 1.
        self.keys = pars.keys()
        self.keys.sort()
        if self.keys not in self._SBkeys:
            import sys
            print 'Not all (or too many) parameters were defined!'
            sys.exit()
        self._baseProfile.__init__(self)
        self.vmap = {}
        self.pars = pars
        for key in self.keys:
            try:
                v = self.pars[key].value
                self.vmap[key] = self.pars[key]
            except:
                self.__setattr__(key,self.pars[key])
        self.setPars()
        self.name = name
        self.convolve = convolve


    def __setattr__(self,key,value):
        if key=='pa':
            self.__dict__['pa'] = value
            if value is not None:
                self.__dict__['theta'] = value*pi/180.
        elif key=='theta':
            if value is not None:
                self.__dict__['pa'] = value*180./pi
            self.__dict__['theta'] = value
        elif key=='logamp':
            if value is not None:
                self.__dict__['amp'] = 10**value
        else:
            self.__dict__[key] = value


    def setPars(self):
        for key in self.vmap:
            self.__setattr__(key,self.vmap[key].value)


class SBTemplate:
    def __init__(self, name, pars, filters, normrange=(4000., 5000.)):
        if 'amp' not in pars.keys() and 'logamp' not in pars.keys():
            pars['amp'] = 1.
        self.keys = pars.keys()
        self.keys.sort()
        if self.keys not in self._SBkeys:
            import sys
            print 'Not all (or too many) parameters were defined!'
            sys.exit()
        self._baseProfile.__init__(self)
        self.vmap = {}
        self.pars = pars
        for key in self.keys:
            try:
                v = self.pars[key].value
                self.vmap[key] = self.pars[key]
            except:
                self.__setattr__(key,self.pars[key])
        self.setPars()
        self.name = name

        splinedic = {}
        rangedic = {}

        wavmax = 0.
        for band in filters:
            splinedic[band] = splrep(filters[band][0], filters[band][1])
            rangedic[band] = (filters[band][0][0], filters[band][0][-1])
            if rangedic[band][1] > wavmax:
                wavmax = rangedic[band][1]

        self.filtsplines = splinedic
        self.filtranges = rangedic

        f = open(tempdir+'CWWSB4.list', 'r')
        lines = f.readlines()
        f.close()

        templates = []
        
        for line in lines:
            line = line.rstrip()
            f = open(tempdir+line, 'r')
            tab = np.loadtxt(f)
            f.close()

            wav = tab[:, 0]
            flambda = tab[:, 1]

            wrange = wav < wavmax

            norm = 1./np.median(flambda[(wav > normrange[0]) & (wav < normrange[1])])

            templates.append((wav[wrange], norm * flambda[wrange]))

        self.templates = templates

    def __setattr__(self,key,value):
        if key=='pa':
            self.__dict__['pa'] = value
            if value is not None:
                self.__dict__['theta'] = value*pi/180.
        elif key=='theta':
            if value is not None:
                self.__dict__['pa'] = value*180./pi
            self.__dict__['theta'] = value
        elif key=='logamp':
            if value is not None:
                self.__dict__['amp'] = 10**value
        else:
            self.__dict__[key] = value


    def setPars(self):
        for key in self.vmap:
            self.__setattr__(key,self.vmap[key].value)

class Sersic(SBModel,SBProfiles.Sersic):
    _baseProfile = SBProfiles.Sersic
    _SBkeys = [['amp','n','pa','q','re','x','y'],
                ['logamp','n','pa','q','re','x','y'],
                ['amp','n','q','re','theta','x','y'],
                ['logamp','n','q','re','theta','x','y']]

    def __init__(self,name,pars,convolve=0):
        SBModel.__init__(self,name,pars,convolve)

    def getMag(self,amp,zp):
        from scipy.special import gamma
        from math import exp,pi
        n = self.n
        re = self.re
        k = 2.*n-1./3+4./(405.*n)+46/(25515.*n**2)
        cnts = (re**2)*amp*exp(k)*n*(k**(-2*n))*gamma(2*n)*2*pi
        return cnts2mag(cnts,zp)

    def Mag(self,zp):
        return self.getMag(self.amp,zp)


class SersicTemplate(SBTemplate, SBProfiles.SersicTemplate):
    _baseProfile = SBProfiles.SersicTemplate
    _SBkeys = [['amp','n','pa','q','re', 'tn', 'x','y', 'zs'],
                ['logamp','n','pa','q','re', 'tn', 'x','y', 'zs'],
                ['amp','n','q','re','theta', 'tn', 'x','y', 'zs'],
                ['logamp','n','q','re','theta', 'tn', 'x','y', 'zs']]

    def __init__(self, name, pars, filters):
        SBTemplate.__init__(self, name, pars, filters)

    def getMag(self,amp,zp):
        from scipy.special import gamma
        from math import exp,pi
        n = self.n
        re = self.re
        k = 2.*n-1./3+4./(405.*n)+46/(25515.*n**2)
        cnts = (re**2)*amp*exp(k)*n*(k**(-2*n))*gamma(2*n)*2*pi
        return cnts2mag(cnts,zp)

    def Mag(self,zp):
        return self.getMag(self.amp,zp)

    def scale(self, band2, band1, madau=True):
        zs = self.zs

        a2 = self.tn%1
        a1 = 1. - a2
        n1 = int(np.floor(self.tn))
        n2 = int(np.ceil(self.tn))

        bands_here = [band2, band1]
        temp_coeffs = [a1, a2]
        temp_ind = [n1, n2]

        temp_fnus = []
        for n in temp_ind:

            flambda = self.templates[n][1].copy()
            wav = self.templates[n][0].copy()

            wobs = wav * (1. + zs)

            fnus_here = []
            for band in bands_here:
                wrange = (wobs > self.filtranges[band][0]) & (wobs < self.filtranges[band][1])
                wband = wobs[wrange]
                fband = flambda[wrange]

                if madau:
                    madau_corr = etau_madau(wband, zs)
                    fband *= madau_corr

                fnu = fband * wband**2
                weights = splev(wband, self.filtsplines[band])
                fnus_here.append((fnu * weights).sum()/weights.sum())

            temp_fnus.append(fnus_here)

        ratio = (a1 * temp_fnus[0][0] + a2 * temp_fnus[1][0]) / (a1 * temp_fnus[0][1] + a2 * temp_fnus[1][1])

        return ratio

class Sersic_wboxyness(SBModel, SBProfiles.Sersic_wboxyness):
    _baseProfile = SBProfiles.Sersic_wboxyness
    _SBkeys = [['amp','b4', 'n','pa','q','re','x','y'],
                ['b4', 'logamp','n','pa','q','re','x','y'],
                ['amp','b4', 'n','q','re','theta','x','y'],
                ['b4', 'logamp','n','q','re','theta','x','y']]

    def __init__(self,name,pars,convolve=0):
        SBModel.__init__(self,name,pars,convolve)

    def getMag(self,amp,zp):
        from scipy.special import gamma
        from math import exp,pi
        n = self.n
        re = self.re
        k = 2.*n-1./3+4./(405.*n)+46/(25515.*n**2)
        cnts = (re**2)*amp*exp(k)*n*(k**(-2*n))*gamma(2*n)*2*pi
        return cnts2mag(cnts,zp)

    def Mag(self,zp):
        return self.getMag(self.amp,zp)


class Spiral(SBModel, SBProfiles.Spiral):
    _baseProfile = SBProfiles.Spiral
    _SBkeys = [['amp', 'bar', 'disk', 'h', 'omega', 'pa', 'q', 'ra', 'x', 'y']]

    def __init__(self, name, pars, convolve=0):
        SBModel.__init__(self, name, pars, convolve)

    def getMag(self, amp, zp):
	cnts = amp
        return cnts2mag(cnts, zp)

    def Mag(self, zp):
        return self.getMag(self.amp, zp)


class Ring(SBModel, SBProfiles.Ring):
    _baseProfile = SBProfiles.Ring
    _SBkeys = [['amp', 'hi', 'ho', 'pa', 'q', 'rr', 'x', 'y']]

    def __init__(self, name, pars, convolve=0):
        SBModel.__init__(self, name, pars, convolve)

    def getMag(self, amp, zp):
	cnts = amp
        return cnts2mag(cnts, zp)

    def Mag(self, zp):
        return self.getMag(self.amp, zp)


class StoneRing(SBModel, SBProfiles.StoneRing):
    _baseProfile = SBProfiles.StoneRing
    _SBkeys = [['amp', 'omega', 'pa', 'q', 'rr', 'smooth', 'spa', 'stone', 'width', 'x', 'y']]

    def __init__(self, name, pars, convolve=0):
        SBModel.__init__(self, name, pars, convolve)

    def getMag(self, amp, zp):
	cnts = amp
        return cnts2mag(cnts, zp)

    def Mag(self, zp):
        return self.getMag(self.amp, zp)


class Arc(SBModel, SBProfiles.Arc):
    _baseProfile = SBProfiles.Arc
    _SBkeys = [['amp', 'hr', 'ht', 'invrc', 'length', 'pa', 'x', 'y']]

    def __init__(self, name, pars, convolve=0):
        SBModel.__init__(self, name, pars, convolve)

    def getMag(self, amp, zp):
        cnts = amp
        return cnts2mag(cnts, zp)

    def Mag(self, zp):
        return self.getMag(self.amp, zp)

class Gauss(SBModel,SBProfiles.Gauss):
    _baseProfile = SBProfiles.Gauss
    _SBkeys = [['amp','pa','q','r0','sigma','x','y']]

    def __init__(self,name,pars,convolve=0):
        if 'r0' not in pars.keys():
            pars['r0'] = None
        SBModel.__init__(self,name,pars,convolve)

    def getMag(self,amp,zp):
        from math import exp,pi
        if self.r0 is None:
            cnts = amp/(2*pi*self.sigma**2)
        else:
            from scipy.special import erf
            r0 = self.r0
            s = self.sigma
            r2pi = (2*pi)**0.5
            cnts = amp*pi*s*(r2pi*r0*(1.+erf(r0/(s*2**0.5)))+2*s*exp(-0.5*r0**2/s**2))
        return cnts2mag(cnts,zp)

    def Mag(self,zp):
        return self.getMag(self.amp,zp)



class BAHBA:
    def __setattr__(self,key,value):
        if key=='logamp':
            if value is not None:
                self.__dict__['amp'] = 10**value
        else:
            self.__dict__[key] = value

    def pixeval(self,xc,yc,dummy1=None,dummy2=None,**kwargs):
        if self.ispix==True:
            return PM.pixeval(self,xc,yc)
        else:
            return GM.pixeval(self,xc,yc)

    def setValues(self):
        self.x = self.values['x']
        self.y = self.values['y']
        if 'amp' in self.keys:
            self.amp = self.values['amp']
        elif self.values['logamp'] is not None:
            self.amp = 10**self.values['logamp']

    def getMag(self,amp,zp):
        return cnts2mag(amp,zp)

    def Mag(self,zp):
        return self.getMag(self.amp,zp)

    def setPars(self,pars):
        for key in self.vmap:
            self.values[self.vmap[key]] = pars[key]
        self.setValues()

