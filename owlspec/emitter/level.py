#!/usr/bin/python

import os
import re
from ..util import parse_spectroscopic_name, eV_in_wn, get_wl
import numpy as np
import pandas as pd

class level():
    def __init__(self, emitter_name, energy, debug=False):
        self.name, self.charge = parse_spectroscopic_name(emitter_name)
        self.spec_name = emitter_name
        from mendeleev import element
        self.emitter = self.particle = element(self.name)
        self.emitter.charge = self.charge
        self.emitter.m = self.emitter.mass
        self.emitter.Ei = self.emitter.ionenergies[1]
        self.charge = self.emitter.charge
        
        from ASDCache import SpectraCache
        nist_cache = SpectraCache()
        nist_levels = nist_cache.levels.fetch(emitter_name, wl_range=(1, 99999))

        # select closest level
        E = nist_levels['Level (cm-1)']/eV_in_wn
        line_idx = (E - energy).abs().idxmin()
        self.E = E.loc[line_idx]

        row = nist_levels.loc[line_idx] 

        self.J = row['J']
        if pd.isna(self.J):
            self.J = None

        try:
            conf_end = row['Configuration'].split('.')
            if "(" in conf_end[-1]:
                conf_end.pop()
            self.l = int(self.l_name_to_num(conf_end[-1][1]))
        except:
            self.l = np.nan
        
        self.g = row['g'] # statistical weight
        if pd.isna(self.g):
            self.g = None

        self.conf = row['Configuration'] # configuration and term string
        
        if 'Lande' in row.keys():
            self.G = row['Lande'] #Lande g factor
        else:
            self.G = None
        
        if pd.isna(self.G):
            self.G = None


    def l_name_to_num(self, name):
        chars = ["s","p","d","f","g","h","i","j"]
        num = chars.index(name)
        return num
    
    
    def get_transitions(self, kind="all", debug=False):
        """Return the transition objects that belong to the energy level
        Type is a string and may be "from", "to" or "all", to select transitions
        from the level, to the level or both."""
        from .transition import transition # needed to avoid circular import

        from ASDCache import SpectraCache
        nist_cache = SpectraCache()
        nist_lines = nist_cache.fetch(self.spec_name, wl_range=(1, 99999))

        transitions_to = []
        transitions_from = []
        transitions_all = []

        upperE = nist_lines['Ek(cm-1)']/eV_in_wn
        lowerE = nist_lines['Ei(cm-1)']/eV_in_wn

        ## upper level
        mask = (upperE - self.E).abs() <= 0.000001
        matches = nist_lines.loc[mask.fillna(False)]

        wl = get_wl(matches)
        for this_wl in wl:
            this_transition = transition(self.spec_name, this_wl)
            transitions_from.append(this_transition)                    
            transitions_all.append(this_transition)

        ## lower level
        mask = (lowerE - self.E).abs() <= 0.000001
        matches = nist_lines.loc[mask.fillna(False)]

        wl = get_wl(matches)
        for this_wl in wl:
            this_transition = transition(self.spec_name, this_wl)
            transitions_to.append(this_transition)                    
            transitions_all.append(this_transition)

        if kind=="all":
            return transitions_all
        if kind=="from":
            return transitions_from        
        if kind=="to":
            return transitions_to
        
        return None


    def get_lifetime(self, debug=False):
        Aiks = []
        for transition in self.get_transitions(kind="from", debug=debug):
            if transition.Aik:
                Aiks.append(transition.Aik)
        sum_Aik = np.sum(Aiks)
        return 1/sum_Aik
