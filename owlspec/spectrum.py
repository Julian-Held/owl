#!/bin/python3
import numpy as np
import pandas as pd
from .util import *
import re
import warnings


class spectrum():
    """ Calculates complete or partial emission spectra for the emitter 
    specified by <emitter name> inr the wavelength range defined by 
    <wl_range>.

    <emitter_name>: Name of the emitting species in spectroscopic notation,
    e.g. "Ar I" for argon neutrals or Fe IV for tripply charged iron.
    <wl_range>: Tuple of lower and upper wavelength of the spectrum."""

    def __init__(self, emitter_name, wl_range=None):
        self.name, self.charge = parse_spectroscopic_name(emitter_name)
        self.spec_name = get_spectroscopic_name(self.name, self.charge)
        from mendeleev import element
        self.emitter = self.particle = element(self.name)
        self.emitter.charge = self.charge
        self.emitter.m = self.emitter.mass
        self.emitter.Ei = self.emitter.ionenergies[1]
        self.charge = self.emitter.charge
        
        self.wl_range = wl_range
        self.linedata = None


    def get_linedata(self):
        if self.wl_range is not None:
            wl_low = self.wl_range[0]
            wl_high = self.wl_range[-1]
        else:
            wl_low = 1
            wl_high = 9999

        try:
            from ASDCache import SpectraCache
            nist_cache = SpectraCache()
            nist_lines = nist_cache.fetch(self.spec_name, wl_range=(wl_low, wl_high))           
        except:
            return []

        return nist_lines


    def get_spectrum(self, x, width=0.02, mu=0.2, norm=False, min_int=-1, min_Aik=-1):
        """ Return simulated spectrum. Lines are pseudo Voigt with
        the set width (FWHM) in nm and form parameter mu (0 = Gauss)  """
        nist_lines = self.get_linedata()
        spectrum = self.table_to_spec_rel(x, nist_lines, width, mu, min_int, min_Aik)
        if norm == True and np.min(spectrum) > 0:
            spectrum = spectrum/np.max(spectrum)
        return spectrum


    def get_LTE_spectrum(self, x, Te, width=0.02, mu=0.2, norm=False, min_int=-1, min_Aik=-1):
        """ Return simulated spectrum with LTE line intnsities in units 
        proportional (!) to Photons/s (NOT W/cm²s)
        Lines are pseudo Voigt with the set width (FWHM) in nm and form 
        parameter mu (0 = Gauss).
        Set norm=True to normalize the intensity to 1 for easier fitting.
        """
        nist_lines = self.get_linedata()
        spectrum = self.table_to_spec_LTE(x, nist_lines, Te, width, mu, min_int, min_Aik)
        if norm == True and np.min(spectrum) > 0:
            spectrum = spectrum/np.max(spectrum)
        return spectrum
        
        
    def get_ident_spectrum(self, min_int=-1, min_Aik=-1, norm=False):
        """ Return simulated spectrum. Marks the wavelenght position
        with thin lines. """
        nist_lines = self.get_linedata()
        x, y = self.table_to_ident(nist_lines, min_int, min_Aik)
        if norm == True and len(y) > 0:
            y = y/np.max(y)
        return x,y


    def get_ident_spectrum_LTE(self, Te, min_int=-1, min_Aik=-1, norm=False):
        """ Return simulated spectrum. Marks the wavelenght position
        with thin lines. """   
        nist_lines = self.get_linedata()
        x, y = self.table_to_ident_LTE(nist_lines, Te, min_int, min_Aik)
        if norm == True and len(y) > 0:
            y = y/np.max(y)
        return x,y
        

    def table_to_spec_LTE(self, x, nist_lines, Te, width=0.1, mu=0.5, min_int=-1, min_Aik=-1):
        spectrum = np.zeros(len(x))

        wl = get_wl(nist_lines)
        Aik = nist_lines['Aki(s^-1)']

        if min_int > 0:
            rel_int = nist_lines['intens']
        else:
            rel_int = 1

        valid = (wl.notna()) & (Aik.notna()) & (rel_int > min_int) & (Aik > min_Aik)
        # these_lines = lines[valid]

        Aik = Aik[valid]
        wl = wl[valid]
        gi = nist_lines['g_k'][valid]
        Ei = nist_lines['Ek(cm-1)'][valid]/eV_in_wn

        intensity = Aik*gi*np.exp(-Ei/Te)
        for this_wl, this_int in zip(wl, intensity):
            profile = this_int*psd_voigt_function(x, this_wl, width, mu)
            spectrum = spectrum + profile
    
        return spectrum
    
    
    def table_to_spec_rel(self, x, nist_lines, width=0.02, mu=0, min_int=-1, min_Aik=-1):
        spectrum = np.zeros(len(x))

        wl = get_wl(nist_lines)
        rel_int = nist_lines['intens']

        if min_Aik > 0:
            Aik = nist_lines['Aki(s^-1)']
        else:
            Aik = 1

        valid = (wl.notna()) & (rel_int.notna()) & (rel_int > min_int) & (Aik > min_Aik)

        if min_Aik > 0:
            Aik = Aik[valid]

        wl = wl[valid]
        rel_int = rel_int[valid]

        for this_wl, this_rel_int in zip(wl, rel_int):
            profile = this_rel_int*psd_voigt_function(x, this_wl, width, mu)
            spectrum = spectrum + profile
                
        return spectrum

        
    def table_to_ident(self, nist_lines, min_int=-1, min_Aik=-1):
        x = []
        dx = 1e-6
        y = []

        wl = get_wl(nist_lines)
        rel_int = nist_lines['intens']

        for idx, row in nist_lines.iterrows():
            this_wl = wl.iloc[idx]
            this_rel_int = rel_int.iloc[idx] 
            if min_Aik > 0:
                this_Aik = row['Aki(s^-1)']
            else:
                this_Aik = 1

            if pd.notna(this_rel_int) and this_rel_int > min_int \
                    and pd.notna(this_Aik) and this_Aik > min_Aik and pd.notna(this_wl):

                x.append(this_wl-dx)
                x.append(this_wl)
                x.append(this_wl+dx)
                y.append(0)
                y.append(this_rel_int)
                y.append(0)
            
        return np.array(x),np.array(y)  
    
        
    def table_to_ident_LTE(self, nist_lines, Te, min_int=-1, min_Aik=-1):
        x = []
        dx = 1e-6
        y = []

        wl = get_wl(nist_lines)
        rel_int = nist_lines['intens']

        for idx, row in nist_lines.iterrows():
            this_wl = wl.iloc[idx]
            this_rel_int = rel_int.iloc[idx] 
            this_Aik = row['Aki(s^-1)']
            this_gi = row['g_k']
            this_Ei = row['Ek(cm-1)']/eV_in_wn

            if min_int > 0:
                this_rel_int = row['intens']
            else:
                this_rel_int = 1
                
            intensity = this_Aik*this_gi*np.exp(-this_Ei/Te)

            if pd.notna(this_rel_int) and this_rel_int > min_int \
                    and pd.notna(intensity) and pd.notna(this_wl):

                x.append(this_wl-dx)
                x.append(this_wl)
                x.append(this_wl+dx)
                y.append(0)
                y.append(intensity)
                y.append(0)

        return np.array(x),np.array(y)  
        
        
