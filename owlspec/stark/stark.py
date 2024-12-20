#!/usr/bin/python

import numpy as np
from ..util import interpol, lorentz_function
from .gigosos_loader import gigosos_loader
from .gigosos_he_loader import gigosos_he_loader
from .griem import griem

class stark():
    def __init__(self, transition, ion_pert = None):
        self.transition = transition
        self.fast = False
        self.pert = ion_pert

    def supported(self):
        wl = self.transition.wl
        x = np.linspace(wl-10, wl+10, 10000)
        ne = 2e21
        profile = None
        try:
            profile = self.get_profile(x, ne, Te=30000, pert=None)
        except:
            return False
        if profile is not None:
            return True
        else:
            return False

    def get_profile(self,x, ne, Te=None, pert=None):
        middle_wl = x[int(len(x)/2)]
        ele = self.transition.emitter.symbol

        ################ Hydrogen ########
        if ele == "H":
            if round(self.transition.wl,0) == 656:
                # Gigosos et al, Spectrochimica Acta Part B 58 (2003) 1489–1504
                if self.fast or not Te:
                    w = ((ne/1e23)**(0.67965)) * 1.098
                    return lorentz_function(x,middle_wl,w)
                else:
                    this_loader = gigosos_loader(self.transition)
                    gigosos_x,y = this_loader.load(ne, Te, pert)
                    gigosos_x = gigosos_x + middle_wl # to nm
                    y = interpol(gigosos_x,y,x)
                    return y

            if round(self.transition.wl, 0) == 486:
                # Gigosos et al, Spectrochimica Acta Part B 58 (2003) 1489–1504
                if self.fast or not Te:
                    w = ((ne/1e23)**(0.68116)) * 4.8
                    return lorentz_function(x,middle_wl,w)
                else:
                    this_loader = gigosos_loader(self.transition)
                    gigosos_x,y = this_loader.load(ne, Te, pert)
                    gigosos_x = gigosos_x + middle_wl # to nm
                    y = interpol(gigosos_x,y,x)
                    return y
                
            if round(self.transition.wl, 0) == 434:
                # Gigosos et al, Spectrochimica Acta Part B 58 (2003) 1489–1504
                this_loader = gigosos_loader(self.transition)
                gigosos_x,y = this_loader.load(ne, Te, pert)
                gigosos_x = gigosos_x + middle_wl # to nm
                y = interpol(gigosos_x,y,x)
                return y

        ################ Helium ##########
        if ele == "He":
            if round(self.transition.wl,0) == 447:
                # Gigosos et al, A&A 503, 293–299 (2009)
                # doi: 10.1051/0004-6361/200912243
                this_loader = gigosos_he_loader(self.transition)
                gigosos_x,y = this_loader.load(ne, Te, pert)
                gigosos_x = gigosos_x + middle_wl # to nm
                y = interpol(gigosos_x,y,x)
                return y

            if round(self.transition.wl, 0) == 492:
                # Lara et al, A&A 542, A75 (2012)
                # doi: 10.1051/0004-6361/201219123
                this_loader = gigosos_he_loader(self.transition)
                gigosos_x,y = this_loader.load(ne, Te, pert)
                gigosos_x = gigosos_x + middle_wl # to nm
                y = interpol(gigosos_x,y,x)
                return y

        ################ Others: Griem ##########
        if ele == "O" or ele == "Ar":
            this_griem = griem(self.transition)
            w,d = this_griem.get_width_shift(ne, Te)
            # shift is ignored for now
            if w>0:
                return lorentz_function(x,middle_wl,w)



    def get_width(self, ne, Te=None, pert=None):
        ele = self.transition.element
        if ele == "O" or ele == "Ar":
            this_griem = griem(self.transition)
            w,d = this_griem.get_width_shift(ne, Te)
            return w

    def get_shift(self, ne, Te=None, pert=None):
        ele = self.transition.element
        if ele == "O" or ele == "Ar" :
            this_griem = griem(self.transition)
            w,d = this_griem.get_width_shift(ne, Te)
            return d


