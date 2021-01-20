# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 14:10:58 2020

@author: Jérémie Marlhens
Title: rbf_mmd
"""

################################################################################
### Quadratic-time MMD with Gaussian RBF kernel

import numpy as np


def rbf(X, Y, sigma=0):
    a=np.shape(X)
    b=np.shape(Y)
    

    G=np.square(X).sum(axis=1)
    H=np.square(Y).sum(axis=1)
    
    Q=np.tile(G,(a[0],1)).transpose()
    R=np.tile(H,(b[0],1))
    
    H= Q+R-2*np.dot(X,Y.transpose())
    
    return np.exp(-H/(2*sigma**2))
    
    
def nmmd2_rbf(X,Y,sigma) :
    a=np.shape(X)
    return np.sum(rbf(X, X, sigma)+rbf(Y, Y, sigma)-2*rbf(X, Y, sigma))/a[0]


