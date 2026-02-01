#!/usr/bin/python3

import numpy as np
from scipy import constants as const
from ..util import interpol
import os
import sys
from platformdirs import user_data_dir
import requests
from pathlib import Path
import tarfile
from io import BytesIO
data_folder447 = os.path.join(user_data_dir("owl-OES", "owl-OES"), "Gigosos2009He")
data_folder492 = os.path.join(user_data_dir("owl-OES", "owl-OES"), "Lara2012He")

            
class gigosos_he_loader():
    def __init__(self, transition):
        self.transition = transition


    def download_profiles(redownload = False):
        if os.path.exists(data_folder447+"/table06.txt") and \
           os.path.exists(data_folder492+"/table06.txt") and redownload==False:
            return
        print("Downloading Stark broadening data tables for helium.")
        URL447 = "http://cdsarc.u-strasbg.fr/ftp/J/A+A/503/293/tab.tar.gz"
        URL492 = "http://cdsarc.u-strasbg.fr/ftp/J/A+A/542/A75/tab.tar.gz"
        try:
            response = requests.get(URL447)
            zip_file447 = response.content
            response = requests.get(URL492)
            zip_file492 = response.content
        except Exception as e:
            print(e)
            print("ERROR: could not download helium Stark broadening tables")
            
        try:
            Path(data_folder447).mkdir(parents=True, exist_ok=True)
            Path(data_folder492).mkdir(parents=True, exist_ok=True)
            with tarfile.open(name=None, fileobj=BytesIO(zip_file447)) as zip_ref:
                zip_ref.extractall(data_folder447)
            with tarfile.open(name=None, fileobj=BytesIO(zip_file492)) as zip_ref:
                zip_ref.extractall(data_folder492)
        except Exception as e:
            print(e)
            print("ERROR: could now save helium Stark broadening tables to disk")
            

    def load(self, ne, Te, pert):
        if pert:
           emitter_m = self.transition.particle.m
           reduced_m = (emitter_m * pert.m)/(emitter_m + pert.m)
           reduced_m = reduced_m/const.u
           mu = reduced_m * Te/pert.T
        else:
           mu = 1

        return self.load_stark_profile(ne, mu, Te*const.eV/const.k)


    def closest(self,val,array):
        return array[np.argmin(np.abs(array-val))]


    def slightly_smaller(self,val,array):
        array = array[array < val]
        return array[np.argmin(np.abs(array-val))]


    def slightly_bigger(self,val,array):
        array = array[array > val]
        return array[np.argmin(np.abs(array-val))]


    def mirror_profile(self,profile):
        return np.append(np.array(profile[::-1]), profile)


    def mirror_x(self,x):
        return np.append(np.array(x[::-1])*-1, x)


    def closest_ne(self,ne):
        if round(self.transition.wl,0) == 447:
            ne_stock = 10**np.linspace(21,24,10)
        if round(self.transition.wl,0) == 492:
            ne_stock = 10**np.linspace(20,24,13)
        if ne in ne_stock:
            ne_low = ne_high = ne
        else:
            ne_low  = self.slightly_smaller(ne,ne_stock)
            ne_high = self.slightly_bigger(ne,ne_stock)
        return ne_low,ne_high


    def closest_mu(self,mu):
        if round(self.transition.wl,0) == 447:
            mu_stock = np.array([0.8, 2, 4])
        if round(self.transition.wl,0) == 492:
            mu_stock = np.array([0.8, 2, 4, 10])
        if mu in mu_stock:
            mu_low = mu_high = mu
        elif mu < mu_stock[0]:
            mu_low = mu_high = mu_stock[0]
        elif mu > mu_stock[-1]:
            mu_low = mu_high = mu_stock[-1]
        else:
            mu_low  = self.slightly_smaller(mu,mu_stock)
            mu_high = self.slightly_bigger(mu,mu_stock)
        return mu_low,mu_high


    def get_T_stock(self, ne):
        if round(self.transition.wl,0) == 447:
            ne_stock = 10**np.linspace(21,24,10) # prevent rounding errors
            if ne <= ne_stock[4]:
                T_stock = np.array([5000, 10000, 20000, 40000])
            if ne > ne_stock[4] and ne < ne_stock[8]:
                T_stock = np.array([10000, 20000, 40000])
            if ne >= ne_stock[8]:
                T_stock = np.array([20000, 40000])
                
        if round(self.transition.wl,0) == 492:
            if ne <= 1e22:
                T_stock = np.array([5000, 10000, 20000, 40000])
            if ne > 1e22 and ne <= 1e23:
                T_stock = np.array([10000, 20000, 30000, 40000])
            if ne > 1e23 and ne < 1e24:
                T_stock = np.array([20000, 30000, 40000])
            if ne == 1e24:
                T_stock = np.array([30000, 40000])
                
        return T_stock


    def closest_T(self,T,ne):
        T_stock = self.get_T_stock(ne)
        if T in T_stock:
            T_low = T_high = T
        elif T < T_stock[0]:
            T_low = T_high = T_stock[0]
        elif T > T_stock[-1]:
            T_low = T_high = T_stock[-1]
        else:
            T_low  = self.slightly_smaller(T,T_stock)
            T_high = self.slightly_bigger(T,T_stock)
        return T_low,T_high


    def load_file(self,ne):
        if round(self.transition.wl,0) == 447:
            folder = data_folder447
            filename = "table"
            table_num = str(int(round((np.log10((ne/1e20))*3-1)))).zfill(2)
        if round(self.transition.wl,0) == 492:
            folder = data_folder492
            folder = os.path.join(folder,"He492")
            filename = "table"
            table_num = str(int(round((np.log10((ne/1e20))*3+4)))).zfill(2)

        filename = filename + table_num + ".txt"
        try:
            table = np.loadtxt(os.path.join(folder,filename)).T
        except:
            raise SystemExit('Error: Could not find Stark broadening data tables.')
        return table


    def interpolate_stark_profile(self, x, low_y, high_y, low_val, high_val, val):
        "Creates profile for arbitrary T/mu. Arrays must use same x!"
        if high_val != low_val:
            low_y = np.array(low_y)
            high_y = np.array(high_y)
            mix_factor = ((val - low_val)/(high_val - low_val))#**(1.5)
            return x,(mix_factor * high_y + (1-mix_factor) * low_y)
        else:
            return x,high_y

    def interpolate_stark_profile_ne(self,low_x, high_x, low_y, high_y, low_val, high_val, val):
        "Creates profile for arbitrary ne. X can be interpolated"
        low_y = interpol(low_x, low_y, high_x) # all to high_x
        if high_val != low_val:
            low_y = np.array(low_y)
            high_y = np.array(high_y)
            mix_factor = ((val - low_val)/(high_val - low_val))
            return high_x,(mix_factor * high_y + (1-mix_factor) * low_y)
        else:
            return high_x,high_y


    def load_stark_profile(self,ne_, mu_, T_):
        if round(self.transition.wl,0) == 447:
            mu_stock = np.array([0.8, 2, 4])
        if round(self.transition.wl,0) == 492:
            mu_stock = np.array([0.8, 2, 4, 10])

        x = None
        profiles_ne = []
        vals_ne = []
        for ne in self.closest_ne(ne_):
            table = self.load_file(ne)
            profiles_mu = []
            vals_mu = []
            for mu in self.closest_mu(mu_):
                T_stock = self.get_T_stock(ne)
                T_low,T_high = self.closest_T(T_,ne)
                pos_low = (list(T_stock).index(T_low)+1) + (len(T_stock) * list(mu_stock).index(mu))
                pos_high = (list(T_stock).index(T_high)+1) + (len(T_stock) * list(mu_stock).index(mu))
                _,y_low = table[0], table[pos_low]
                high_x,y_high = table[0], table[pos_high]
                x,y = self.interpolate_stark_profile(high_x,y_low,y_high,T_low,T_high,T_)
                profiles_mu.append((x,y))
                vals_mu.append(mu)

            x,y = self.interpolate_stark_profile(
                profiles_mu[1][0],profiles_mu[0][1],
                profiles_mu[1][1],vals_mu[0],vals_mu[1],mu_)

            profiles_ne.append((x,y))
            vals_ne.append(ne)

        x,y = self.interpolate_stark_profile_ne(
            profiles_ne[0][0],profiles_ne[1][0],profiles_ne[0][1],
            profiles_ne[1][1],vals_ne[0],vals_ne[1],ne_)

        return x,y



