[![DOI](https://zenodo.org/badge/267619393.svg)](https://zenodo.org/doi/10.5281/zenodo.11002438)

# owl spectroscopic library
  
Library to calculate optical emission spectra from atoms as well as line shapes.

## Installation
Owl is on [pypi](https://pypi.org/project/owlspec/) and can simply be installed as any python package:
```
pip install owlspec
```

For manual installation, you can also simply drop the `owlspec` folder in the folder where your script is located. However, you will still need to install the dependencies (numpy, scipy, mendeleev, astroquery, roman).

### Downloading Stark profiles
Before using the library for Stark broadening of the hydrogen Balmer series as well as helium 447/492 nm, you need to download the precalculated tables by Gigosos et al/Lara et al. 
In case of the hydrogen Balmer series, the neccesary files are hosted by Elsevier only for private non-commercial use (see the [copyright notice](https://github.com/Julian-Held/owl/blob/master/README.md#copyright-notice-for-stark-broadening-tables) from the [publication](https://doi.org/10.1016/S0584-8547(03)00097-1)).

Downloading the data tables for hydrogen Balmer Stark broadeing:
``` python
import owlspec as owl
owl.gigosos_loader.download_profiles()
```

For the helium lines at 447 nm and 492 nm:

``` python
import owlspec as owl
owl.gigosos_he_loader.download_profiles()
```

You can simply place these lines at the top of any script if you like: the files will not be re-downloaded if they are already there. If the files get corrupted and you want to redownload them, you can call the function with `redownload = True`.

## Capabilities
The library has three core capabilities:
1. Fetch information about transitions and levels from NIST.
2. Calculating the broadening of emission lines for a selection of cases.
3. Calculating complete PLTE spectra.

### Levels and transitions data
Owl can be used to access information from the NIST atomic spectra database[^NIST], which we access using [astroquery](https://astroquery.readthedocs.io/en/latest/). That does mean _owl requires internet access to function_. To access information about a transition, create a transition object by supplying the emitter name in spectroscopic notation (e.g. "Ar I" for argon neutrals or "Fe IV" for triply charged iron ions) and the wavelength in nm:
```python
transition = owl.emitter.transition("Ar I", 751.5)
```
The object then contains the information as attributes. For example, the energy of the upper level of the transition can be accessed as `transition.upperE`, the Einstein coefficient for emission as `transition.Aik` and so on. Possible attributes are: `name`, `charge`, `spec_name`, `wl`, `upperE`, `lowerE`, `upperl`, `lowerl`,`Aik`,`upperg`, `lowerg`. Additionally, the levels corresponding to the transition can be accessed as `transition.upper` and `transition.lower`. The levels contain `J` and the Lande g factor `G`.

Alternatively, the energy levels can also be accessed directly by specifying the energy (in eV):
```python
level = owl.emitter.level("Ar I", 11.623)
```
The level contains attributes `E`, `J`, `l` and `G` and `conf`. 

By itself, the level and transition data is not very useful, but would be helpful in, for example, assembling a collisional radiative model. However, the level class contains a method that might be of immediate use: `get_lifetime()` which calculates the radiative lifetime of the level from all transitions listed in NIST. The result is returned in units of seconds.

### Line broadening calculations

Owl supports a range of different line broadening mechanisms that are automatically activated when providing information about the physical situation surrounding the emitters. For example, specifying an electron density will automatically switch on Stark broadening calculations, if they are available for the emitter. We currently support Doppler, Stark, Zeeman, van der Waals and instrumental broadening. Self-resonance broadening is not yet supported. Multiple broadening mechanisms are combined by numerical convolution of the individual profiles.

#### Stark broadening
Since Stark broadening is different for each transition and each emitter species and no generalized theory is available, owl only supports a few selected transitions.

Owl uses the Stark broadening calculations of Gigosos et al for H alpha, beta and gamma, as well as He 447.1 nm and He 492.2 nm. The precalculated tables are downloaded from the publishers when the respective functions are first executed. The space in between the calculated datapoints is interpolated as advised by the respective authors.

Additionally, we have support for O 777 nm, Ar 810.369 nm and Ar 738.398 nm based on Griems tabulated constants[^1]. Extending the library with more data from Griems calculations for other transitions would be trivial, but doing everything would be a tremendous amount of busywork. Thus, if you need support for any specific line, let me know and I would be happy to put it in.

| Transition | Wavelength | Source                                                                                                                                         |
| ---------- | ---------- | ---------------------------------------------------------------------------------------------------------------------------------------------- |
| Hα         | 656.3 nm   | [Gigosos et al 2003 _Spectrochimica Acta Part B: Atomic Spectroscopy_ **58** 1489–504](https://doi.org/10.1016/S0584-8547(03)00097-1)          |
| Hβ         | 486.1 nm   | [Gigosos et al 2003 _Spectrochimica Acta Part B: Atomic Spectroscopy_ **58** 1489–504](https://doi.org/10.1016/S0584-8547(03)00097-1)          |
| Hγ         | 434.0 nm   | [Gigosos et al 2003 _Spectrochimica Acta Part B: Atomic Spectroscopy_ **58** 1489–504](https://doi.org/10.1016/S0584-8547(03)00097-1)          |
| He I       | 447.1 nm   | [Gigosos et al 2009 _A&A_ **503** 293–9](https://doi.org/10.1051/0004-6361/200912243)                                                          |
| He I       | 492.2 nm   | [Lara N et al 2012  _A&A_ **542** A75](https://dx.doi.org/10.1051/0004-6361/201219123)                                                         |
| Ar I       | 810.369 nm | [H.R. Griem: Spectral Line Broadening by Plasmas](https://shop.elsevier.com/books/spectral-line-broadening-by-plasmas/griem/978-0-12-302850-1) |
| Ar I       | 738.398 nm | [H.R. Griem: Spectral Line Broadening by Plasmas](https://shop.elsevier.com/books/spectral-line-broadening-by-plasmas/griem/978-0-12-302850-1) |

[^1]: H.R. Griem: Spectral Line Broadening by Plasmas

#### Zeeman splitting
Splitting due to magnetic fields is calculated analytically following the books of Cowan[^Cowan] as well as Condon and Shortley[^CS]. The necessary base data (as the Lande g factors) are obtained from the NIST atomic spectra database [^NIST]. Transitions for which the base data is missing from NIST are not currently supported. Polarization of the different components is calculated, but is currently not exposed to the user. Let me know if you need that feature.

[^Cowan]: Cowan R D 1981 _The theory of atomic structure and spectra_ (Berkeley: University of California Press)
[^CS]: Condon E U and Shortley G H 1979 _The theory of atomic spectra_ (Cambridge: Univ. Pr)
[^NIST]: https://www.nist.gov/pml/atomic-spectra-database

#### Doppler broadening
Doppler broadening currently only supports Maxwell VDFs.

#### van der Waals broadening
Van der Waals broadening is calculated using the equations from Konjević [^Konj]. The polarizability of atoms surrounding the emitters is automatically obtained using [mendeleev](https://mendeleev.readthedocs.io/en/stable/). 

[^Konj]: Konjević  1999 _Physics Reports_ **316** 339–401

#### Instrumental broadening
Instrumental broadening is included either using a pseudo Voigt function with user-specified width and shape parameter or by passing a function that takes the x-axis and the central wavelength position as arguments. The example for Zeeman broadening below shows how to do that.

### Spectra simulations

The library queries the NIST atomic spectra database[^NIST] to automatically obtain base data about transitions and levels. This can be used to calculate complete spectra. For example, to obtain the spectrum of chromium neutrals between 300 nm and 500 nm in partial local thermal equilibrium at 3 eV, all you need to do is:

```python
spec = owl.spectrum('Cr I', wl_range=[300,500])
y = spec.get_LTE_spectrum(x, Te=3, width=0.2, mu=0.5, norm=True) 
```

This capabilities can be used to quickly identify unknown lines in measurements or fit measured spectra to obtain an excitation temperature in cases where Boltzmann plots are difficult due to insufficient resolution. Please note that the calculated spectra are in units proportional to photons/second, _not_ W/(sr cm²).


## Examples

The following examples show most of what owl is currently capable of. All examples should be completely copy-pasteable, if you have the library and matplotlib installed.

### H beta Stark broadening

```python
#!/bin/python3
import numpy as np
import matplotlib.pyplot as plt
import owlspec as owl
owl.gigosos_loader.download_profiles() # Download data tables if neccesary

cw = central_wavelength = 486
x = np.linspace(cw-4, cw+6.5, 3000)

transition1 = owl.emitter.transition("H I", cw) # Hydrogen atom emission, H beta
# H I emitters are 1000 K hot, and sourrounded by argon neutrals and ions 
# and 30000 K hot electrons
line1 = owl.emission_line(transition1, cw, pert="Ar I", T=1000, Te=30000)
# setting T switches on Doppler broadening
# Te and the perturber will be used for Stark broadening below

plt.figure()
y1 = line1.get_profile(x, ne=5e21, ng=1e20) # ng for van der Waals broadening..
y2 = line1.get_profile(x, ne=2e22, ng=1e20) # ..can be set here or in emission_line
y3 = line1.get_profile(x, ne=5e22, ng=1e20)

plt.plot(x,y1/np.max(y1), label=r"n$_e$ = $5 \times 10^{21}$ m$^{-3}$")
plt.plot(x,y2/np.max(y2), label=r"n$_e$ = $2 \times 10^{22}$ m$^{-3}$")
plt.plot(x,y3/np.max(y3), label=r"n$_e$ = $5 \times 10^{22}$ m$^{-3}$")
plt.legend()
plt.xlabel("wavelength / nm")
plt.ylabel("normalized intensity")
```    
<img src="https://github.com/user-attachments/assets/438e9602-8020-46d0-afc1-0b69697638ed" width=400px>


### van der Waals and Doppler broadening of a helium line

```python
#!/bin/python3
import numpy as np
import matplotlib.pyplot as plt
import owlspec as owl

plt.figure() # van der Waals broadening
cw = central_wavelength = 501.567
transition1 = owl.emitter.transition("He I", cw) # Helium atom emission
x = np.linspace(cw-0.1, cw+0.15, 3000)

# He I emitters are 500 K hot, and sourrounded by argon neutrals
line = owl.emission_line(transition1, cw, pert="Ar I", T=500)

y1 = line.get_profile(x, ng=1e24) # ng is the neutral density of the perturber (Ar)
y2 = line.get_profile(x, ng=5e24) # setting ng switches on vdW broadening
y3 = line.get_profile(x, ng=1e25) # T (defined as 500 above) siwtches on Doppler
plt.plot(x,y1/np.max(y1), label=r"n$_g$ = $1 \times 10^{24}$ m$^{-3}$")
plt.plot(x,y2/np.max(y2), label=r"n$_g$ = $5 \times 10^{24}$ m$^{-3}$")
plt.plot(x,y3/np.max(y3), label=r"n$_g$ = $1 \times 10^{25}$ m$^{-3}$")
plt.legend()
plt.xlabel("wavelength / nm")
plt.ylabel("normalized intensity")


plt.figure() # Doppler broadening
y1 = line.get_profile(x, T=1000) # Setting T switches on Doppler broadening
y2 = line.get_profile(x, T=10000)

plt.plot(x,y1/np.max(y1), label="T$_g$ = 1000 K")
plt.plot(x,y2/np.max(y2), label="T$_g$ = 10000 K")
plt.legend()
plt.xlabel("wavelength / nm")
plt.ylabel("normalized intensity")
```

<img src="https://github.com/mimurrayy/owl/assets/3911345/0aab3190-e863-4644-aab7-266676493d01" width=400px><img src="https://github.com/user-attachments/assets/918085e6-7256-4122-9084-357030170310" width=400px>


### Zeeman splitting of Ar I 706.722 nm

```python
#!/bin/python3
import numpy as np
import matplotlib.pyplot as plt
import owlspec as owl

def instr(x, xc): # Function to use for instrumental broadening
    return owl.util.psd_voigt(x, xc, w=0.003, mu=0.2)

plt.figure() # Show Zeeman splitting
cw = central_wavelength = 706.722
x = np.linspace(cw-0.05, cw+0.05, 3000)
transition = owl.emitter.transition("Ar I", cw)
line = owl.emission_line(transition, cw, instr_func=instr)


y1 = line.get_profile(x, B = 0.2)
y2 = line.get_profile(x, B = 0.8)

plt.plot(x,y1/np.max(y1), label="B = 0.2 T")
plt.plot(x,y2/np.max(y2), label="B = 0.8 T")
plt.legend()
plt.xlabel("wavelength / nm")
plt.ylabel("normalized intensity")
```

<img src="https://github.com/user-attachments/assets/066727ec-ef88-4e58-9ba2-b61e35e0eec0" width=400px>

### Using pyplas plasma objects

 Here we use [pyplas](https://github.com/mimurrayy/pyplas) plasma objects which allow us to describe a more complex physical situation, i.e. where species, densities and temperatures differ between emitting species, the background gas (affecting vdW) and the ions (affecting Stark).

```python
#!/bin/python3
import numpy as np
import matplotlib.pyplot as plt
import owlspec as owl
import pyplas
owl.gigosos_loader.download_profiles()

n0 = 5e21 # Electron/Ion density
electrons = pyplas.electrons(n0, T=35000) # 3 eV electrons
ions = pyplas.ions("Ar", n0, T=1000) # 1000 K argon ions 

# Background gas is He, can set any two of n, T and P
neutrals = pyplas.neutrals("He", P=800, T=500) 
plasma = pyplas.plasma(electrons, ions, neutrals) # create plasma object holding all of the information

cw = central_wavelength = 486
transition = owl.emitter.transition("H I", cw)

# T is here now only the emitter temperature, perturber temps are set by plasma
line1 = owl.emission_line(transition, cw, plasma=plasma, T=1000)

x = np.linspace(cw-4, cw+6.5, 3000)
y1 = line1.get_profile(x)
plasma.ne = 2e22 # the plasma object stays connected to the line object.
y2 = line1.get_profile(x)
plasma.ne = 5e22
y3 = line1.get_profile(x)

plt.figure()
plt.plot(x,y1/np.max(y1), label=r"n$_e$ = $5 \times 10^{21}$ m$^{-3}$")
plt.plot(x,y2/np.max(y2), label=r"n$_e$ = $2 \times 10^{22}$ m$^{-3}$")
plt.plot(x,y3/np.max(y3), label=r"n$_e$ = $5 \times 10^{22}$ m$^{-3}$")
plt.legend()
plt.xlabel("wavelength / nm")
plt.ylabel("normalized intensity")
```

<img src="https://github.com/user-attachments/assets/850c566b-e504-4124-9a97-d54aedee3286" width=400px>


### Inspecting the different line broadening contributions
After the calculation of an individual emission line profile is perfromed, the different components are saved to a 'profiles' dictonary in the emission_line object.
You can access these profiles to ensure that the individual calculations are perfromed correctly and to judge the contribution of the seperate mechanisms.

```python
#!/bin/python3
import numpy as np
import matplotlib.pyplot as plt
import owlspec as owl
import pyplas
owl.gigosos_loader.download_profiles()

n0 = 5e20 # Electron/Ion density
electrons = pyplas.electrons(n0, T=35000) # 3 eV electrons
ions = pyplas.ions("Ar", n0, T=1000) # 1000 K argon ions 
neutrals = pyplas.neutrals("Xe", P=30000, T=300) 
plasma = pyplas.plasma(electrons, ions, neutrals)

cw = central_wavelength = 486

transition = owl.emitter.transition("H I", cw)
line1 = owl.emission_line(transition, cw, plasma=plasma, T=5000)

x = np.linspace(cw-0.2, cw+0.2, 3000)
complete = line1.get_profile(x) # run the calculation ...
Stark = line1.profiles['Stark'] # ... now the different components are saved in the 'profiles' dict in the line object
fine_x =  line1.profiles['x'] # internally, a finer x grid is used, so you need that for plotting
Doppler = line1.profiles['Doppler']
vdW = line1.profiles['vdW']
# not shown here: Instrumental profile accessed with line1.profiles['instrument']

plt.figure()
plt.plot(x,complete/np.max(complete), label='All')
plt.plot(fine_x,Stark/np.max(Stark), '--', label='Stark')
plt.plot(fine_x,Doppler/np.max(Doppler), '--', label='Doppler')
plt.plot(fine_x,vdW/np.max(vdW), '--', label='van der Waals')
plt.legend()
plt.xlabel("wavelength / nm")
plt.ylabel("normalized intensity")
```

<img src="https://github.com/user-attachments/assets/2de7bdca-f7ad-4cbb-8716-690aa5c21a4c" width=400px>

### Switching broadening mechanisms on/off
Instead of relying on automatic decisions made by the library, you can also decide to explicitly switching broadening mechanisms on or off.
The advantage is that owlspec will then notify you if you forgot to supply information which would have otherwise just silently disabled a mechanism.

```python
#!/bin/python3
import numpy as np
import matplotlib.pyplot as plt
import owlspec as owl
import pyplas
owl.gigosos_loader.download_profiles()

n0 = 5e20 # Electron/Ion density
electrons = pyplas.electrons(n0, T=35000) # 3 eV electrons
ions = pyplas.ions("Ar", n0, T=1000) # 1000 K argon ions 
neutrals = pyplas.neutrals("Xe", P=30000, T=300) 
plasma = pyplas.plasma(electrons, ions, neutrals)

cw = central_wavelength = 486

transition = owl.emitter.transition("H I", cw)
line1 = owl.emission_line(transition, cw, plasma=plasma, T=5000)

x = np.linspace(cw-0.2, cw+0.2, 3000)
complete = line1.get_profile(x, Stark_on=True, Doppler_on=True, Instr_on=True, vdW_on=True)
vdW = line1.get_profile(x, Stark_on=False, Doppler_on=False, Instr_on=False, vdW_on=True)
Stark = line1.get_profile(x, Stark_on=True, Doppler_on=False, Instr_on=False, vdW_on=False)
Doppler = line1.get_profile(x, Stark_on=False, Doppler_on=True, Instr_on=False, vdW_on=False)


plt.figure()
plt.plot(x,complete/np.max(complete), label='All')
plt.plot(x,Stark/np.max(Stark), '--', label='Stark')
plt.plot(x,Doppler/np.max(Doppler), '--', label='Doppler')
plt.plot(x,vdW/np.max(vdW), '--', label='van der Waals')
plt.legend()
plt.xlabel("wavelength / nm")
plt.ylabel("normalized intensity")
```

<img src="https://github.com/user-attachments/assets/2de7bdca-f7ad-4cbb-8716-690aa5c21a4c" width=400px>


### Stick spectrum for the identification of emission lines
In order to identify unknown lines in a measured spectrum, it's useful to compare the measurement to lines reported in NIST using their reported relative intensities. To this end, `owl.spectrum.get_ident_spectrum()` can be used, which marks the line position and intensity with a single dot, which is much faster than calculating realistic line shapes. Intensities according to PLTE are also supported, by using `owl.spectrum.get_ident_spectrum_LTE(Te)`, instead.

```python
#!/bin/python3
import numpy as np
import matplotlib.pyplot as plt
import owlspec as owl

plt.figure() # Ident spectrum

 # Select all Ti lines between 600 and 900 nm listed by NIST
spec = owl.spectrum('Ti I', wl_range=[310,415])
spec2 = owl.spectrum('Ar II', wl_range=[310,415])

# simulate spectra with a given line braodening (w, mu) and normalized to the maximum
# min_int and min_Aik allow filtering for NIST-reported relative intensity and Aik value
x,y = spec.get_ident_spectrum(min_int=200, min_Aik=1e6) 
x2,y2 = spec2.get_ident_spectrum() 
plt.plot(x,y, '-', label="Ti I")
plt.plot(x2,y2, '-', label="Ar II")
plt.legend()
plt.xlabel("wavelength / nm")
plt.ylabel("relative intensity (NIST)")
```

<img src="https://github.com/mimurrayy/owl/assets/3911345/c57b40c2-d596-40c9-bf70-101847030393" width=400px>


### Simulation of a complete PLTE spectrum

The example below shows how to calculate a complete LTE spectrum. The function is fast enough that it can be used to fit measured data to obtain an excitation temperature.

```python
#!/bin/python3
import numpy as np
import matplotlib.pyplot as plt
import owlspec as owl

plt.figure() # PLTE Spectrum simulation
x = np.linspace(310,415,10000)

 # Select all Ti lines between 600 and 900 nm listed by NIST
spec = owl.spectrum('Ti I', wl_range=[310,415])

# simulate spectra with a given line braodening (w, mu) and normalized to the maximum
y = spec.get_LTE_spectrum(x, Te=3, width=0.2, mu=0.5, norm=True) 
y2 = spec.get_LTE_spectrum(x, Te=1, width=0.2, mu=0.5, norm=True)
plt.plot(x,y, '-', label="T$_e$ = 3 eV")
plt.plot(x,y2, '-', label="T$_e$ = 1 eV")
plt.legend()
plt.xlabel("wavelength / nm")
plt.ylabel("normalized intensity")
```

<img src="https://github.com/mimurrayy/owl/assets/3911345/99354512-8a80-411b-93be-2e4d89b42ce9" width=400px>



## How to cite
The library has a DOI provided by Zenodo and can be cited similar to: 
- Julian Held: (2024) owl spectroscopic library (v0.3.0) [https://doi.org/10.5281/zenodo.11002438](https://doi.org/10.5281/zenodo.11002438)

The library relies on data from the NIST atomic spectra database, which can be cited as:
- Kramida, A., Ralchenko, Yu., Reader, J. and [NIST ASD Team](https://physics.nist.gov/PhysRefData/ASD/index.html#Team) (2023). _NIST Atomic Spectra Database_ (version 5.11). Available: [https://physics.nist.gov/asd](https://physics.nist.gov/asd) [Sat Apr 20 2024]. National Institute of Standards and Technology, Gaithersburg, MD. DOI: [https://doi.org/10.18434/T4W30F](https://doi.org/10.18434/T4W30F)

Up to date citation information for NIST ASD is provided here: https://physics.nist.gov/PhysRefData/ASD/Html/verhist.shtml

When using owl for line broadening calculations, please make sure to cite the sources of the underlying data listed in this readme. 


## Copyright notice for Stark broadening tables
The Stark broadening tables are distributed by Elsevier und the following copyright notice contained within the publication [Gigosos et al 2003 _Spectrochimica Acta Part B: Atomic Spectroscopy_ **58** 1489–504](https://doi.org/10.1016/S0584-8547(03)00097-1):

_The data files are copyrighted by the authors(s). Readers of Spectrochimica Acta Electronica are permitted by the publisher Elsevier B.V., to copy the materialfor their own private, non-commercial use, and to run the programs according to the instructions provided by the authors. No charge for any copies may be requested, neither may the program or any modified version of it be sold or used for commercialpurposes. Those who wish to use the data files for commercial purposes should contact the corresponding author at the address given on the hardcopy paper._



