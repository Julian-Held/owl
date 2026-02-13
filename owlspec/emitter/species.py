#!/usr/bin/python

import numpy as np
import pyplas
from scipy import constants as const


class species(pyplas.species):
    def __init__(self, element, T=None, n=None, P=None, q=0):
        if element:
            super().__init__(element, T=T, n=n, P=P, q=q)
            self.element_info(element)
            self.ideal_gas(n, T, P)
            self.charge = int(round(self._q/const.e + 1, 0))
            self.Ei = self.ionization_energy/const.eV
            self.symbol = self.name = self.element.symbol
            self.dipole_polarizability = self.element.dipole_polarizability

    @property
    def q(self):
        return self._q

    @q.setter
    def q(self, new_value):
        self._q = new_value
        self.element_info(self.element.symbol)


class perturber(species):
    def __init__(self, element, T=None, n=None, P=None, q=0):
        super().__init__(element, T=T, n=n, P=P, q=q)