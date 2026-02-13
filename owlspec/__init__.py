#!/usr/bin/python
from . import emitter
from . import util
from .emission_line import *
from .spectrum import *
from .stark import gigosos_loader, gigosos_he_loader

def download_stark_tables(redownload=False):
    gigosos_loader.download_profiles(redownload)
    gigosos_he_loader.download_profiles(redownload)
