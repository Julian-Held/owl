#!/usr/bin/python3

import numpy as np
from ..util import *
from .base_spectrometer import *

class triax(base_spectrometer):
    def __init__(self, cw, order=1, grating_constant=50, camera="princeton"):
        self.cw = cw
        self.order = order
        self.grating_constant = grating_constant
        self.camera = camera.lower()
        self.x = self.make_x_scale(cw,order,1024)
        self.fine_x = self.make_x_scale(cw,order,64*1024)


    def make_x_scale(self, cw=None, order=1, size=977):
        """Calculates the wavelength scale. Can use bigger sizes to
        create a better resolution for the plots of the results.
        Only works for the 50 line/cm grating for now."""

        if self.camera == "princeton":
            c0 = 5.78
            c1 = 0.688
            x = np.linspace(0,size-1,size)
            x = c0 + cw + (-size/2 + x) * c1
            return x

        elif self.camera == "andor":
            c0 = -2.5
            c1 = 0.345
            x = np.linspace(0,size-1,size)
            x = c0 + cw + (-size/2 + x) * c1
            return x


    def instrument_function(self, x, xc, transition=None, params = None):
        if params:
            w,mu = params
        elif transition:
           w,mu = self.get_instrumental_params(transition)
        y = psd_voigt(x,xc,w,mu)
        return y


    def get_instrumental_params(self, transition):
        "Database for instrumental functions for lines."
        if self.order == 1 and self.grating_constant == 50:
            return (2.48216407e+00, 3.63312977e-01)


    # alias
    instr = instrument_function
    instrumental_function = instrument_function
    instrumental_profile = instrument_function