#!/bin/python3
import numpy as np
from scipy import constants as const
from .util import *
import os
from .emitter import transition
import mendeleev 
from astroquery.nist import Nist
import astropy.units as u
import re
import warnings


class spectrum():
    """ Holds information about a complete emission spectrum of the given
    particles in the wavelength range defined by <wl_range> or the
    owl.spectrometer object.

    particles: list of owl.emitter.particle objects
    spectrometer: owl.spectrometer object
    wl_range: tuple or list of lower and upper wavelength of the spectrum"""

    def __init__(self, emitter_name, spectrometer=None, wl_range=None):
        self.name, self.charge = parse_spectroscopic_name(emitter_name)
        self.spec_name = get_spectroscopic_name(self.name, self.charge)
        self.emitter = self.particle = mendeleev.element(self.name)
        self.emitter.charge = self.charge
        self.emitter.m = self.emitter.mass
        self.emitter.Ei = self.emitter.ionenergies[1]
        self.charge = self.emitter.charge
        
        self.spectrometer = spectrometer
        self.wl_range = wl_range
        self.linedata = None
        if spectrometer and not wl_range:
            if spectrometer.x is not None:
                self.wl_range = (spectrometer.x[0],spectrometer.x[-1])



    def get_linedata(self):
        if self.wl_range is not None:
            wl_low = self.wl_range[0]
            wl_high = self.wl_range[-1]
        else:
            wl_low = 0
            wl_high = 9999
            
        nist_lines = Nist.query(wl_low*u.nm, wl_high*u.nm, 
                                linename=self.spec_name, wavelength_type='vac+air')

        return nist_lines


    def get_spectrum(self,x=None, width=0.02, mu=0.2, min_int=-1, min_Aik=-1):
        """ Return simulated spectrum. Lines are pseudo Voigt with
        the set width (FWHM) in nm and form parameter mu (0 = Gauss)  """
        try:
            if not x:
                x = self.spectrometer.x
        except:
            x = x
            
        nist_lines = self.get_linedata()
        spectrum = self.table_to_spec_rel(x, nist_lines, width, mu, min_int, min_Aik)

        return spectrum


    def get_LTE_spectrum(self, x, Te, width=0.02, mu=0.2, norm=False, min_int=-1, min_Aik=-1):
        """ Return simulated spectrum with LTE line intnsities in units 
        proportional (!) to Photons/s (NOT W/cm²s)
        Lines are pseudo Voigt with the set width (FWHM) in nm and form 
        parameter mu (0 = Gauss).
        Set norm=True to normalize the intensity to 1 for easier fitting.
        """
        try:
            if not x:
                x = self.spectrometer.x
        except:
            x = x
        nist_lines = self.get_linedata()
        spectrum = self.table_to_spec_LTE(x, nist_lines, Te, width, mu, min_int, min_Aik)
        if norm == True:
            spectrum = spectrum/np.max(spectrum)
        return spectrum
        
        
    def get_ident_spectrum(self, min_int=-1, min_Aik=-1):
        """ Return simulated spectrum. Marks the wavelenght position
        with thin lines. """
            
        nist_lines = self.get_linedata()
        x, y = self.table_to_ident(nist_lines, min_int, min_Aik)
        
        return x,y


    def get_ident_spectrum_LTE(self, Te, min_int=-1, min_Aik=-1):
        """ Return simulated spectrum. Marks the wavelenght position
        with thin lines. """
            
        nist_lines = self.get_linedata()
        x, y = self.table_to_ident_LTE(nist_lines, Te, min_int, min_Aik)
        
        return x,y
        

    def table_to_spec_LTE(self, x, nist_lines, Te, width=0.1, mu=0.5, min_int=-1, min_Aik=-1):
        spectrum = np.zeros(len(x))
        for line in  nist_lines:   
            c1 = str(line['Observed'])
            # astroquery does not filter out headings in the middle of the table
            if 'Observed' in c1 or "Wavelength" in c1 or "nm" in c1:
                continue
            
            wl =  line['Observed']
            Aik = line['Aki']
            
            if not np.ma.is_masked(wl) and not np.ma.is_masked(Aik):
                
                wl =  str(wl)
                Aik = str(Aik)
                wl = re.sub(r'[^\d.]+', '', wl) # removes non-number chars
                # left/right strips non-number chars (preserves the 'e' in 5e7 
                Aik = re.sub(r'(^[^\d.]+)|([^\d.]+$)', '', Aik) 
                
                if min_int > 0:
                    with warnings.catch_warnings(): # ignore masked element warning
                        warnings.simplefilter("ignore", category=UserWarning)
                        rel_int = line['Rel.']
                        rel_int = re.sub(r'[^\d.]+', '', rel_int)
                        rel_int =  float(rel_int)
                else:
                    rel_int = 1
    
                try: # conversion to float can still fail...
                    wl = float(wl)
                    Aik = float(Aik)
                    gi = line['gi   gk'].split('-')[1]
                    Ei = line['Ei           Ek'].split('-')[1]
                    Ei = float(Ei)
                    gi = float(gi)
    
                    if rel_int > min_int and Aik > min_Aik:
                        intensity = Aik*gi*np.exp(-Ei/Te)
                        profile = intensity*psd_voigt_function(x, wl, width, mu)
                        spectrum = spectrum + profile
                except Exception as err:
                    print(err)
                    pass
    
        return spectrum
    
    
    def table_to_spec_rel(self, x, nist_lines, width=0.02, mu=0, min_int=-1, min_Aik=-1):
        spectrum = np.zeros(len(x))
        for line in  nist_lines:   
            rel_int = line['Rel.']
            wl =  line['Observed']
            if not np.ma.is_masked(rel_int) and not np.ma.is_masked(wl):
                rel_int = str(rel_int) # is byte type initally
                wl =  str(wl)
                rel_int = re.sub(r'[^\d.]+', '', rel_int)
                wl = re.sub(r'[^\d.]+', '', wl) # removes non-number chars
                try: # conversion to float can still fail...
                    rel_int = float(rel_int)
                    wl = float(wl)
                    if min_Aik > 0:
                        with warnings.catch_warnings(): # ignore masked element warning
                            warnings.simplefilter("ignore", category=UserWarning)
                            Aik =  float(line['Aki'])
                    else:
                        Aik = 1
                        
                    if rel_int > min_int and Aik > min_Aik:
                        profile = rel_int*psd_voigt_function(x, wl, width, mu)
                        spectrum = spectrum + profile
                except:
                    pass
                
        return spectrum

        
    def table_to_ident(self, nist_lines, min_int=-1, min_Aik=-1):
        x = []
        dx = 1e-6
        y = []
        for line in  nist_lines:   
            rel_int = line['Rel.']
            wl =  line['Observed']
            if not np.ma.is_masked(rel_int) and not np.ma.is_masked(wl):
                rel_int = str(rel_int) # is byte type initally
                wl =  str(wl)
                rel_int = re.sub(r'[^\d.]+', '', rel_int)
                wl = re.sub(r'[^\d.]+', '', wl) # removes non-number chars
                try: # conversion to float can still fail...
                    rel_int = float(rel_int)
                    wl = float(wl)
                    if min_Aik > 0:
                        with warnings.catch_warnings(): # ignore masked element warning
                            warnings.simplefilter("ignore", category=UserWarning)
                            Aik =  float(line['Aki'])
                    else:
                        Aik = 1
                        
                    if rel_int > min_int and Aik > min_Aik:
                        x.append(wl-dx)
                        x.append(wl)
                        x.append(wl+dx)
                        y.append(0)
                        y.append(rel_int)
                        y.append(0)
                        
                except:
                    pass
            
        return np.array(x),np.array(y)  
    
        
    def table_to_ident_LTE(self, nist_lines, Te, min_int=-1, min_Aik=-1):
        x = []
        dx = 1e-6
        y = []
        for line in  nist_lines:   
            c1 = str(line['Observed'])
            # astroquery does not filter out headings in the middle of the table
            if 'Observed' in c1 or "Wavelength" in c1 or "nm" in c1:
                continue
            
            wl =  line['Observed']
            Aik = line['Aki']
            
            if not np.ma.is_masked(wl) and not np.ma.is_masked(Aik):
                
                wl =  str(wl)
                Aik = str(Aik)
                wl = re.sub(r'[^\d.]+', '', wl) # removes non-number chars
                # left/right strips non-number chars (preserves the 'e' in 5e7 
                Aik = re.sub(r'(^[^\d.]+)|([^\d.]+$)', '', Aik) 
                
                if min_int > 0:
                    with warnings.catch_warnings(): # ignore masked element warning
                        warnings.simplefilter("ignore", category=UserWarning)
                        rel_int =  float(line['Rel.'])
                else:
                    rel_int = 1

                try: # conversion to float can still fail...
                    wl = float(wl)
                    Aik = float(Aik)
                    gi = line['gi   gk'].split('-')[1]
                    Ei = line['Ei           Ek'].split('-')[1]
                    Ei = float(Ei)
                    gi = float(gi)
                    if rel_int > min_int and Aik > min_Aik:
                        intensity = Aik*gi*np.exp(-Ei/Te)
                        x.append(wl-dx)
                        x.append(wl)
                        x.append(wl+dx)
                        y.append(0)
                        y.append(intensity)
                        y.append(0)
                except:
                    pass

        return np.array(x),np.array(y)  
        
        
