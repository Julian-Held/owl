#!/bin/python3
import numpy as np
from scipy.signal import fftconvolve as convolve
from scipy import constants as const
from scipy import interpolate
from .util import zeeman, parse_spectroscopic_name, \
                        doppler_maxwell, gauss_function, psd_voigt
from . import stark
from . import vdW
from .emitter.species import species, perturber
from .util import parse_spectroscopic_name, get_spectroscopic_name
from copy import copy


def _parse_pert(pert, T=None):
    """Convert a spectroscopic name to a species; pass other objects through."""
    if isinstance(pert, str):
        name, charge = parse_spectroscopic_name(pert)
        pert = species(name, charge*const.e)
        if T:
            pert.T = T
    return pert


def _enabled(flag, available):
    """Resolve a True/False/'auto' switch given whether 'auto' may enable it."""
    return bool(available) if flag == 'auto' else bool(flag)


class emission_line():
    def __init__(self, transition, wl, plasma=None, instr_func = None, w = None, mu = 0.0,
                 T = None, B = None, ne = None, Te = None, ng = None, pert = None, 
                 Instr_on='auto', Doppler_on='auto', Stark_on='auto', Zeeman_on='auto', 
                 vdW_on='auto'):
        
        """  Class to calculate the shape of an emission line based on the 
        physical situation. Providing information on densities and temperatures
        automatically switches on different broadening mechanisms.

        Supported broadening mechanisms: Doppler, Stark, Zeeman, van der Waals
        and instrumental broadening. Stark broadening is only supported for
        some transitions: H-alpha to gamma, He 447.1 nm, He 492.2 nm, O 777 nm,
        Ar 810.369 nm and Ar 738.398 nm.

        <transition>: An owl.emitter.transition object for the specific line.
        <wl>: Center wavelength in nm. Does not need to be the theoretical wl
              of the transition. Use to move simulation relative to measurement.
              Defaults to transition wavelength if not specified.
        <plasma>: pyplas plasma object describing the physical situation. Can
                  be used instead of seperate parameters (T, ne, Te, ng, ...).
                  Needs electron and ion densities and temperatures for Stark.
                  Needs neutral density and temperatures for vdW.
                  Needs magnetic field strength for Zeeman.
                  Not used for Doppler/instrumental broadening.
        <instr_func>: Instrumental function of the spectrometer. Must be a 
                      function that take parameters x and xc (wavelength axis 
                      and central position) and return the profile as a 
                      numpy array.
        <w>: Width of the instrumental function in nm. Only if instr_func=None.
        <mu>: Shape parameter of the instrumental fucntion. 0=Gauss, 1=Lorentz.
        Use together with w or use instr_func instead. Only if instr_func=None.
        <T>: Emitter temperature in Kelvin. Used for vdW and Doppler broadening.
             If both plasma and T are specified, we assume T to be the emitter
             temperature and plasma.Tg=plasma.n.T to be the perturber 
             temperature (otherwise, they are assumed to be the same). 
        <B>: Magnetic field strength for Zeeman splitting. Units of Tesla.
             Use if plasma = None. Superseeds plasma values if both are used.
        <ne>: Electron density for Stark broadening. Units of m^-3. 
              Use if plasma = None. Superseeds plasma values if both are used.
        <Te>: Electron temperature for stark broadening. Units of Kelvin (!).
              Has only minor impact on the results and can usually be guessed.
              Use if plasma = none. Superseeds plasma values if both are used.
        <ng>: Neutral density for van der Waals broadening. Units of m^-3.
             This is the density of species sourrounding the emitter.
        <pert>: Perturber species for vdW and Stark boradening. 
                String in spectroscopic notation, e.g. 'O I' or 'Ar II' or
                owlspec.emitter.species object. This is the ion or atom species 
                surrounding the emitter particles. Note that you cannot specify 
                the perturber temperature when using a string, so we use the 
                emitter temperature set with T. 
                Use if plasma = None.
        <Instr_on>, <Doppler_on>, <Stark_on>, <Zeeman_on>, <vdW_on>: 
        One of True, False 'auto'.
        Turns claculation of respective broadening mechanism on or off. Defaults
        to 'auto' which automatically enables calculation if possible.
        """

        self.transition = transition
        self.wl = wl
        if not wl:
            self.wl = transition.wavelength
        self.plasma = plasma
        self.instr, self.w, self.mu = instr_func, w, mu
        particle = transition.particle
        self.m  = particle.m
        self.Ei = particle.Ei
        self.T, self.B, self.ne, self.Te, self.ng  = T, B, ne, Te, ng
        self.pert = _parse_pert(pert, T)
        self.Instr_on, self.Doppler_on = Instr_on, Doppler_on
        self.Zeeman_on, self.vdW_on, self.Stark_on = Zeeman_on, vdW_on, Stark_on

        self.profiles = dict.fromkeys(['Doppler', 'Zeeman', 'Stark', 'vdW', 'instrument'])

        if isinstance(pert, (species, perturber)) and T and not pert.T:
            pert.T = T


    def get_profile(self, x, A = 1,  wl=None, plasma=None, instr_func = None, 
                    w = None, mu = None, T = None, B = None, ne = None, 
                    Te = None, ng = None, pert = None,
                 Instr_on=None, Doppler_on=None, Stark_on=None, Zeeman_on=None, 
                 vdW_on=None):

        """  Return emission line profile for the physical situation speciefied
        here or in the emission_line object. Values set here superseed object
        variables.

        <x>: numpy array with the wavelength axis in nm.
        <A>: Line intensity (= line area). Defaults to 1.
        <wl>: Center wavelength in nm. Does not need to be the theoretical wl
              of the transition. Use to move simulation relative to measurement.
              Defaults to transition wavelength if not specified.
        <plasma>: pyplas plasma object describing the physical situation. Can
                  be used instead of seperate parameters (T, ne, Te, ng, ...).
                  Needs electron and ion densities and temperatures for Stark.
                  Needs neutral density and temperatures for vdW.
                  Needs magnetic field strength for Zeeman.
                  Not used for Doppler/instrumental broadening.
        <instr_func>: Instrumental function of the spectrometer. Must be a 
                      function that take parameters x and xc (wavelength axis 
                      and central position) and return the profile as a 
                      numpy array.
        <w>: Width of the instrumental function in nm. Only if instr_func=None.
        <mu>: Shape parameter of the instrumental fucntion. 0=Gauss, 1=Lorentz.
        Use together with w or use instr_func instead. Only if instr_func=None.
        <T>: Emitter temperature in Kelvin. Used for vdW and Doppler broadening.
             If both plasma and T are specified, we assume T to be the emitter
             temperature and plasma.Tg=plasma.n.T to be the perturber 
             temperature (otherwise, they are assumed to be the same). 
        <B>: Magnetic field strength for Zeeman splitting. Units of Tesla.
             Use if plasma = None. Superseeds plasma values if both are used.
        <ne>: Electron density for Stark broadening. Units of m^-3. 
              Use if plasma = None. Superseeds plasma values if both are used.
        <Te>: Electron temperature for Stark broadening. Units of Kelvin (!).
              Has only minor impact on the results and can usually be guessed.
              Use if plasma = none. Superseeds plasma values if both are used.
        <ng>: Neutral density for van der Waals broadening. Units of m^-3.
             This is the density of species sourrounding the emitter.
        <pert>: Perturber species for vdW and Stark boradening. 
                String in spectroscopic notation, e.g. 'O I' or 'Ar II' or
                owlspec.emitter.species object. This is the ion or atom species 
                surrounding the emitter particles. Note that you cannot specify 
                the perturber temperature when using a string, so we use the 
                emitter temperature set with T. 
                Use if plasma = None.
        <Instr_on>, <Doppler_on>, <Stark_on>, <Zeeman_on>, <vdW_on>: 
        One of True, False 'auto'.
        Turns claculation of respective broadening mechanism on or off. Defaults
        to 'auto' which automatically enables calculation if possible.
        """

        def pick(value, default):
            return default if value is None else value

        wl = pick(wl, self.wl)
        plasma = pick(plasma, self.plasma)
        instr_func = pick(instr_func, self.instr)
        w = pick(w, self.w)
        mu = pick(mu, self.mu)
        T = pick(T, self.T)
        B = pick(B, self.B)
        ne = pick(ne, self.ne)
        Te = pick(Te, self.Te)
        ng = pick(ng, self.ng)
        Instr_on = pick(Instr_on, self.Instr_on)
        Doppler_on = pick(Doppler_on, self.Doppler_on)
        Stark_on = pick(Stark_on, self.Stark_on)
        Zeeman_on = pick(Zeeman_on, self.Zeeman_on)
        vdW_on = pick(vdW_on, self.vdW_on)

        # self.pert is already parsed; only a per-call override needs parsing
        pert = self.pert if pert is None else _parse_pert(pert, T)

        vdW_pert = None
        Stark_pert = None

        if not pert and ng and T:
            vdW_pert = self.transition.particle
            vdW_pert.T = T
            vdW_pert.n = ng
        
        if not pert and ne and T:
            Stark_pert = self.transition.particle
            Stark_pert.T = T
            Stark_pert.n = ne

        if plasma:
            if plasma.ions:
                Stark_pert = copy(plasma.ions)
            if plasma.neutrals:
                vdW_pert = copy(plasma.neutrals)
        elif pert:
            Stark_pert = copy(pert)
            vdW_pert = copy(pert)
            
        if T and ng:
            if T and not vdW_pert.T:
                vdW_pert.T = T
            if ng:
                vdW_pert.n = ng

        if T and ne:
            if T and not Stark_pert.T:
                Stark_pert.T = T
            if ne and not Stark_pert.n:
                Stark_pert.n = ne

        if plasma:
            if plasma.ne and not ne:
                ne = plasma.ne
            if plasma.Te and not Te:
                Te = plasma.Te
            
            if plasma.Tg and not T:
                if self.transition.emitter.q > 0:
                    T = plasma.ions.T
                else:
                    T = plasma.Tg

            if not B:
                B = plasma.B
        

        self.profiles = dict.fromkeys(['x', 'y', 'Doppler', 'Zeeman', 'Stark', 'vdW', 'instrument'])
        orig_x = np.copy(x)
        s = 0 # lineshift in nm
        fine_x = interpolate.make_interp_spline(np.linspace(0,1,len(x)), x)
        x = fine_x(np.linspace(0,1,len(x)*50)) # upsample x axis
        self.profiles['x'] = x

        resolution = abs((x[-1]-x[0])/(len(x)-1))
        middle_wl = x[int(len(x)/2)-1]
        components = []

        if _enabled(Doppler_on, T):
            maxwell = doppler_maxwell(x, middle_wl, T, self.m)
            self.profiles['Doppler'] = maxwell/np.max(maxwell)
            components.append(maxwell)

        if _enabled(Stark_on, ne):
            this_stark = stark.stark(self.transition)
            stark_profile = this_stark.get_profile(x, ne, Te, Stark_pert)
            self.profiles['Stark'] = stark_profile/np.max(stark_profile)
            components.append(stark_profile)

        if _enabled(Zeeman_on, B):
            if not (B):
                print('Missing magnetic field strength for Zeeman splitting. Skipping.')
            else:
                t = self.transition
                if t.upperJ and t.lowerJ and t.upperG and t.lowerG:
                    zeeman_pattern = zeeman(x, middle_wl, B, t.upperJ, t.lowerJ, t.upperG, t.lowerG)/resolution
                    self.profiles['Zeeman'] = zeeman_pattern/np.max(zeeman_pattern)
                    components.append(zeeman_pattern)
                else:
                    print('Missing quantum numbers (J) or Lande G values for Zeeman splitting. Skipping.')

        if _enabled(vdW_on, vdW_pert):
            if not (vdW_pert.T and vdW_pert.n):
                if vdW_on==True:
                    print('Missing perturber/neutral density and temperature for van der Waals broadening. Skipping.')
            else:
                this_vdW = vdW.vdW(self.transition, pert=vdW_pert)
                vdW_profile = this_vdW.get_profile(x, vdW_pert.n, vdW_pert.T)
                self.profiles['vdW'] = vdW_profile/np.max(vdW_profile)
                components.append(vdW_profile)
                # s = s + this_vdW.get_shift(x, vdW_pert.n, vdW_pert.T)

        # We let the instrumental profile determine center position.
        # ALL other components are shifted to the middle!
        instr_enabled = _enabled(Instr_on, True)
        if instr_func and instr_enabled:
            instrumental_profile = instr_func(x, wl+s)
        elif w and instr_enabled:
            instrumental_profile = psd_voigt(x, wl+s, w, mu)
        else:
            instrumental_profile = gauss_function(x, wl+s, x[1]-x[0])
        
        self.profiles['instrument'] = instrumental_profile/np.max(instrumental_profile)

        profile = instrumental_profile
        for component in components:
            profile = convolve(profile, component, mode='same') * resolution

        y = A*profile
        lowres_y = np.interp(orig_x, x, y)
        self.profiles['y'] = A*y/np.sum(lowres_y)/resolution
        
        return A*lowres_y/np.sum(lowres_y)/resolution# normalize intensity to A
