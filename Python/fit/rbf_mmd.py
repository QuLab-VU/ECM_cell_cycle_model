# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 14:10:58 2020

@author: Jérémie Marlhens
Title: rbf_mmd
"""

################################################################################
### Quadratic-time MMD with Gaussian RBF kernel

import numpy as np

def rbf(X, Y, sigma):
    Q = np.transpose(np.tile(X.T[:, :, None], (1, 1, Y.shape[0])), (1, 0, 2))
    R = np.transpose( np.tile(Y.T[:, :, None], (1, 1, X.shape[0])), (2, 0, 1))
    return np.exp(-np.sum((Q-R)**2, axis=1) / (2 * sigma ** 2))

def mmd_rbf(X, Y, sigma):
    m = X.shape[0]
    n = Y.shape[0]
    
    Kxx = rbf(X, X, sigma)
    Kyy = rbf(Y, Y, sigma)
    Kxy = rbf(X, Y, sigma)
    
    np.fill_diagonal(Kxx, 0)
    np.fill_diagonal(Kyy, 0)
    
    return (np.sum(Kxx) / (m * (m - 1)) + np.sum(Kyy) / (n * (n - 1)) - 2 * np.sum(Kxy) / (m * n))



if __name__ == "__main__":
    sigma = 1
    X = np.array([np.random.rand(10)]).transpose()
    Y = np.array([np.random.rand(20)]).transpose()
    mmd_rbf(X, Y, sigma)
    
    test_10=0
    for i in range(1000):
        X = np.array([np.random.rand(10)]).transpose()
        Y = 4 * np.array([np.random.rand(10)]).transpose()
        test_10 += mmd_rbf(X, Y, sigma)
    test_10/1000
        
    test_100=0
    for i in range(1000):
        X = np.array([np.random.rand(100)]).transpose()
        Y = 4 * np.array([np.random.rand(100)]).transpose()
        test_100 += mmd_rbf(X, Y, sigma)
    test_100/1000
    
