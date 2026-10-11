#!/bin/python3
import numpy as np
from scipy.signal import fftconvolve as convolve
from scipy import constants as const
from scipy import interpolate
from scipy.signal import filtfilt
from numpy import fft
import roman

__all__ = [
    "interpol",
    "doppler_maxwell",
    "gauss_function","gauss",
    "m_j",
    "zeeman",
    "get_spectroscopic_name",
    "parse_spectroscopic_name",
    "running_mean",
    "downsample",
    "fft_clean",
    "lorentz_function","lorentz",
    "psd_voigt_function","psd_voigt",
    "fft_smooth",
    "deconv", "deconvolution",
    "get_wl",
    "eV_in_wn",
    ]

eV_in_wn = 8065.544  # express energy in eV

def interpol(new_x, y, old_x):
    f = interpolate.interp1d(new_x, y, fill_value=0.0, bounds_error=False)
    return f(old_x)

def doppler_maxwell(x, xc, T, m):
    w = xc/const.c * np.sqrt(8 * const.k * T * np.log(2) / m)
    return gauss_function(x, xc, w)

def gauss_function(x, xc, w):
    "A is Area, w is FWHM, H = y0 + A/w / np.sqrt(0.25 * np.pi / np.log(2))"
    return (1/(w*np.sqrt(np.pi/(4*np.log(2)))) *
        np.exp(-4*np.log(2)*(x-xc)**2/w**2)) # gaussian function copied from origin, w is FWHM

def m_j(j):
    if j == round(j,0):
        return np.append(-np.arange(j+1)[::-1][:-1],np.arange(j+1))
    if j == round(j,0) + 0.5 or  j == round(j,0) - 0.5:
        return np.append(-np.arange(j+0.5)[::-1][:-1],np.arange(j+1.5))-0.5

def zeeman(x, cwl, B, upperJ, lowerJ, upperG, lowerG, side=False):
    bm =  const.physical_constants['Bohr magneton'][0]
    z_peaks = []
    y = np.zeros(len(x))
    J_u = upperJ
    J_l = lowerJ
    counter = 0 # for intensity scaling
    for mu in m_j(J_u):
        for ml in m_j(J_l):
            if abs(mu-ml) <= 1: # transition optical allowed?
                counter = counter + 1
                dE = bm * B * mu * upperG - bm * B * ml * lowerG
                E0 = const.h * const.c / (cwl * 1e-9)
                wl = cwl * (dE/E0 + 1)

                # intensity acording to Condon,Shortley Page 387
                if J_u == J_l: # J -> J transition
                    if mu == ml: # pi transition with delta M = 0
                        intensity = mu**2
                        if side:
                            intensity = 0
                    if mu == (ml + 1): # M -> M - 1
                        intensity = 0.25 * (J_u + mu) * (J_u - mu + 1)
                    if mu == (ml - 1): # M -> M + 1
                        intensity = 0.25 * (J_u - mu) * (J_u + mu + 1)
                if J_u == (J_l - 1): # J -> J + 1 transition
                    if mu == ml:
                        intensity = (J_u + 1)**2 - mu**2
                        if side:
                            intensity = 0
                    if mu == (ml + 1): # M -> M - 1
                        intensity = 0.25 * (J_u - mu + 1) * (J_u - mu + 2)
                    if mu == (ml - 1): # M -> M + 1
                        intensity = 0.25 * (J_u + mu + 1) * (J_u + mu + 2)
                if J_u == (J_l + 1): # J -> J -1 transition
                    if mu == ml:
                        intensity = J_u**2 - mu**2
                        if side:
                            intensity = 0
                    if mu == (ml + 1): # M -> M - 1
                        intensity = 0.25 * (J_u + mu) * (J_u + mu - 1)
                    if mu == (ml - 1): # M -> M + 1
                        intensity = 0.25 * (J_u - mu) * (J_u - mu - 1)
                index = np.argmin(np.abs(x-wl)) # index of line component
                y[index] = y[index] + intensity # in case two transitions hit the same pixel
    y = y/counter
    return y

def get_spectroscopic_name(element, charge):
    name = element.title() # ti -> Ti
    designation = roman.toRoman(round(charge+1))
    return name + " " + designation

def parse_spectroscopic_name(name):
    name = name.strip()
    ele = name.split(' ')[0]
    num = name.split(' ')[-1]
    charge = roman.fromRoman(num) - 1
    return ele.title(), charge

def running_mean(x, N):
    cumsum = np.cumsum(np.insert(x, 0, 0))
    return (cumsum[N:] - cumsum[:-N]) / float(N)

def downsample(y, N):
    end =  N * int(len(y)/N)
    return np.mean(y[:end].reshape(-1, N), 1)

def fft_clean(y, cutoff):
    Hn = np.fft.rfft(y)
    Hn[cutoff:] = 0
    return np.fft.irfft(Hn)

def fft_smooth(y, smoothness):
    n = smoothness # the larger n is, the smoother curve will be
    b = [1.0 / n] * n
    a = 1.0
    return filtfilt(b,a,y)

def lorentz_function(x, xc, w):
    "A is area, w is FWHM, H = y0 + 2*A / (np.pi * w)"
    return ((2*1/np.pi)*(w/(4*(x-xc)**2 + w**2)))

def psd_voigt_function(x, xc, w, mu):
    return (1 * ( mu * (2/np.pi) * (w / (4*(x-xc)**2 + w**2)) +
        (1 - mu) * (np.sqrt(4*np.log(2)) / (np.sqrt(np.pi) * w)) *
        np.exp(-(4*np.log(2)/w**2)*(x-xc)**2) )) # pseudo voidt function copied from origin

def deconv(signal, instr, noise_level):
    """ Remove instrumental profile (or other profile) from measured signal.
    Wiener deconvolution found somewhere on the internet.
    Was checked to be correct by forward calculation using convolution.
    Noise level is 1/SNR with SNR being the signal to noise ratio.
    In practice, increase noise level until result starts looking smushed.
    """
    H = fft.fft(instr)
    lamb=max(instr)*noise_level
    deconvolved = fft.ifftshift(fft.ifft(fft.fft(signal)*
                                         np.conj(H)/(H*np.conj(H) + lamb**2)))
    return np.roll(np.real(deconvolved), -1)

def get_wl(nist_lines):
    """ Get wavelength: try observed is possible, otherwise use Ritz values

        Under 200 nm and above 2000 nm, use vacuum wavelengths instead of air wavelengths """

    wl_vac = nist_lines['obs_wl_vac(nm)']
    wl_vac[wl_vac.isna()] = nist_lines[wl_vac.isna()]['ritz_wl_vac(nm)']

    wl = nist_lines['obs_wl_air(nm)']
    wl[wl.isna()] = nist_lines[wl.isna()]['ritz_wl_air(nm)']

    wl[wl_vac < 200] = wl_vac[wl_vac < 200]
    wl[wl_vac > 2000] = wl_vac[wl_vac > 2000]

    return wl

# alias
lorentz = lorentz_function
psd_voigt = psd_voigt_function
gauss = gauss_function
deconvolution = deconv
