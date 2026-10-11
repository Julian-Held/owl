#!/usr/bin/python

from operator import index
import pandas as pd

import numpy as np
import re
from ..util import parse_spectroscopic_name, get_spectroscopic_name, get_wl, eV_in_wn
from .level import level
from .species import species
from scipy import constants as const



class transition():
    def __init__(self, emitter, wavelength, debug=False):
        """Return the closest transition for the specified emitter to the
        specified wavelength.
        <emitter>: String in spectroscopic notation, e.g. O I or Ar II or
                   owlspec.emitter.species object.
        <wavelength>: Wavelength specifiying the transition in nanometers. 
                      Chooses the closest match from NIST, might need to 
                      specify to picometer precision.
        
        Returns a transition object with entries: upperE, lowerE, upperl, 
        lowerl, upperg, lowerg, Aik, wl and the upper 
        and lower levels of the transition: upper_level/lower_level."""
        
        if isinstance(emitter, str):
            self.name, self.charge = parse_spectroscopic_name(emitter)
            self.spec_name = get_spectroscopic_name(self.name, self.charge)
            self.emitter = self.particle = species(self.name, q=self.charge*const.e)
        
        elif isinstance(emitter, species):
            self.spec_name = get_spectroscopic_name(emitter.name, emitter.charge)
            self.charge = emitter.charge
            self.emitter = self.particle = emitter
               
        # Load data tables from NIST (+- 1 nm around requested wl)
        from ASDCache import SpectraCache
        nist_cache = SpectraCache()
        nist_lines = nist_cache.fetch(self.spec_name, wl_range=(wavelength-1, wavelength +1))  
        
        # select closest line
        wl = get_wl(nist_lines)
        line_idx = (wl - wavelength).abs().idxmin()
        self.wl = wl.loc[line_idx]

        row = nist_lines.loc[line_idx] 

        self.upperE = row["Ek(cm-1)"] / eV_in_wn
        self.lowerE = row["Ei(cm-1)"] / eV_in_wn

        
        upper_conf_end = row['conf_k'].split('.')
        if "(" in upper_conf_end[-1]: # not sure if still needed
            upper_conf_end.pop()
        lower_conf_end = row['conf_i'].split('.')
        if "(" in lower_conf_end[-1]:
            lower_conf_end.pop()
        try:
            self.upperl = self.l_name_to_num(upper_conf_end[-1][1])
            self.lowerl = self.l_name_to_num(lower_conf_end[-1][1]) 
        except:
            self.upperl = None
            self.lowerl = None
        
        self.Aik    = row['Aki(s^-1)']
        self.upperg = float(row['g_k'])
        self.lowerg = float(row['g_i'])

        if pd.isna(self.upperg):
            self.upperg = None
        if pd.isna(self.lowerg):
            self.lowerg = None
        

        self.upper_level, self.lower_level = self.levels()
        self.upper, self.lower = self.upper_level, self.lower_level

    def l_name_to_num(self, name):
        chars = ["s","p","d","f","g","h","i","j"]
        num = chars.index(name)
        return num

    def levels(self):
        upper_level = level(self.spec_name, self.upperE)
        lower_level = level(self.spec_name, self.lowerE)
        self.upperJ = upper_level.J
        self.upperG = upper_level.G
        self.lowerJ = lower_level.J
        self.lowerG = lower_level.G
        return upper_level, lower_level

    