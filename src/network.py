#!/usr/bin/env python

# Load required modules
import numpy as np
import scipy as sp

##########################################################################
# DIFFUSION MEASURES
##########################################################################

def regularized_laplacian(L, lam):
    return np.linalg.inv(np.eye( *np.shape(L) ) + (lam * L))

def rkhs_factor(D):
    e, v = sp.linalg.eigh(D)
    return v.dot(np.diag(np.sqrt(e)))

