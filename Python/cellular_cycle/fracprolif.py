# -*- coding: utf-8 -*-
"""
Created on Sat Aug 22 11:24:52 2020

@author: Grim.M
"""

path="/home/jeremie"



import numpy as np
from scipy.optimize import minimize, fmin_l_bfgs_b, brentq, approx_fprime, least_squares
from scipy.special import erf
from scipy.stats import norm
import scipy.integrate as integrate
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from scipy.optimize import leastsq

#############################################################################

def demg(x, mu=0, sigma=1, lambd=1):
    return (lambd/2) * np.exp((lambd/2)*(2*mu -2*x + lambd*(sigma*sigma))) * (1 - erf((-x+mu+lambd*(sigma**2))/(np.sqrt(2)*sigma)))

def pemg(q, mu = 0, sigma = 1, lambd = 1, lower_tail = True):
    u = (q - mu) * lambd
    sigma1 = lambd * sigma
    p = norm.cdf(u, 0, sigma1) - np.exp(-u + (sigma1 * sigma1)/2 +  np.log(norm.cdf(u, sigma1 * sigma1, sigma1)))
    if lower_tail==False:
        p=1-p
    return p


def qsurvival_nllik(x0):
    mu, sigma, lambd, Q = x0   
    global mitotic_lifespans, censored_lifespans
    return - sum([np.log(1-Q) + np.log(demg(i, mu = mu, sigma = sigma, lambd = lambd)) for i in  mitotic_lifespans]) - sum([np.log(Q+(1-Q)*pemg(i,  mu = mu, sigma = sigma, lambd = lambd, lower_tail = False)) for i in censored_lifespans])
    



def mle(mitotic_lifespans, censored_lifespans):
    x0 = np.array([np.median(mitotic_lifespans), np.std(mitotic_lifespans), 1/np.mean(mitotic_lifespans), 0.5])
    
    def qsurvival_nllik(x0):
        mu, sigma, lambd, Q = x0   
        #global mitotic_lifespans, censored_lifespans
        return - sum([np.log(1-Q) + np.log(demg(i, mu = mu, sigma = sigma, lambd = lambd)) for i in  mitotic_lifespans]) - sum([np.log(Q+(1-Q)*pemg(i,  mu = mu, sigma = sigma, lambd = lambd, lower_tail = False)) for i in censored_lifespans])
        return - sum([np.log(1-Q) + np.log(demg(i, mu = mu, sigma = sigma, lambd = lambd)) for i in  mitotic_lifespans])
    
    bounds = [(1, 100), (0.01, (max(mitotic_lifespans) - min(mitotic_lifespans))), (0.001, 100), (1e-41, 0.999999999999)]
    #
    #m = minimize(qsurvival_nllik, x0, method='nelder-mead', bounds=bounds) #Change the method to minimize
    #x: array([18.89482881,  1.67120341,  0.09735528,  0.59049594])
    #x0: array([25.5       , 10.35847627,  0.03450854,  0.5       ])
    
    #m = minimize(qsurvival_nllik, x0, method="L-BFGS-B", bounds=bounds)
    
    #m = minimize(qsurvival_nllik, x0, method="TNC", bounds=bounds)
    #print(mitotic_lifespans)
    #print(censored_lifespans)
    #print(x0)
    #print(bounds)
    try:
        m = minimize(qsurvival_nllik, x0, method="SLSQP", bounds=bounds)
        return m
    except:
        pass
#        np.save(path + "/fit_data/pydream/frac_prolif_quiescence_DMSO_only_30_parameters/mitotic_lifespans.npy", mitotic_lifespans)
#        np.save(path + "/fit_data/pydream/frac_prolif_quiescence_DMSO_only_30_parameters/censored_lifespans.npy", censored_lifespans)
#        np.save(path + "/fit_data/pydream/frac_prolif_quiescence_DMSO_only_30_parameters/x0.npy", x0)
#        np.save(path + "/fit_data/pydream/frac_prolif_quiescence_DMSO_only_30_parameters/bounds.npy", bounds)
    #m = minimize(qsurvival_nllik, x0, method="trust-constr", bounds=bounds)        


def q_rates_inner(m_x) : #m_x is 
    def integrand(t, v):
        return np.exp(-v*t)*demg(t, mu=m_x[0], sigma=m_x[1], lambd=m_x[2])
    def integr(v) : 
        return integrate.quad(integrand, 0, np.inf, args=(v))[0] - 0.5
    root = brentq(integr, 1e-04, 1)
    d = root
    q = m_x[3] * root/(1 - m_x[3])
    return d, q


def q_rates(m):
    n = m['nit']
#    def q_rates_inner_d(m_x):
#        return q_rates_inner(m_x)[0]
#    def q_rates_inner_q(m_x):
#        return q_rates_inner(m_x)[1]
#    
    rates=q_rates_inner(m['x']) #
#    eps = m['x']/1000
#    gd = approx_fprime(m['x'], q_rates_inner_d, eps)
#    gq = approx_fprime(m['x'], q_rates_inner_q, eps)
#    
    v={}
    v['d']=rates[0]
    v['q']=rates[1]
    #v["stderr"] = c(sqrt(t(gd) %*% vcov(est)[1:3, 1:3] %*% gd/n), sqrt(t(gq) %*% vcov(est) %*% gq/n))
    #v["n"] = n
    #v["df"] = n - 2
    return v


def prolif_q(t, d, q, ad, aq, x0, q0 = 0, tol = 1e-06) :
    if abs(d - q) < tol :
        if abs(ad - aq) < tol :
            return np.exp(-aq * t) * (q0 + x0 * q * t)
        else:
            return np.exp(-aq * t)/(ad - aq) * (x0 * q * (1 - exp((aq - ad) * t)) + (ad - aq) * q0)
    else :
        beta = d - q + aq - ad
        return np.exp(-aq * t)/beta * (q * (x0 * (np.exp(beta * t) - 1) - q0) + (beta + q) * q0)

    
def prolif_d(t, d, q, ad, aq, x0, q0 = 0, tol = 1e-06):
    return x0 * np.exp((d - q - ad) * t)


def split_tot(t, d1, d2, q1, q2, ad1, aq1, ad2, aq2, x0, q0=0, br = 20) :
    qbr = prolif_q(br, d1, q1, ad1, aq1, x0, q0=q0) #
    dbr = prolif_d(br, d1, q1, ad1, aq1, x0, q0=q0) #n
    
    before = prolif_d(t, d1, q1, ad1, aq1, x0, q0=q0) + prolif_q(t, d1, q1, ad1, aq1, x0, q0=q0)
    
    after = prolif_d(t - br, d2, q2, ad2, aq2, dbr, q0=qbr) + prolif_q(t - br, d2, q2, ad2, aq2, dbr, q0=qbr)
    
    answer = before
    if t > br:
        answer=after
    return answer

def split_q(t, d1, d2, q1, q2, ad1, aq1, ad2, aq2,x0, q0=0, br = 20) :
    qbr = prolif_q(br, d1, q1, ad1, aq1, x0, q0=q0) #
    dbr = prolif_d(br, d1, q1, ad1, aq1, x0, q0=q0) #n
    before = prolif_q(t, d1, q1, ad1, aq1, x0, q0=q0)
    after = prolif_q(t - br, d2, q2, ad2, aq2, dbr, q0=qbr)
    answer = before
    if t > br:
        answer=after
    return answer

def split_d(t, d1, d2, q1, q2, ad1, aq1, ad2, aq2,x0, q0=0, br = 20) :
    qbr = prolif_q(br, d1, q1, ad1, aq1, x0, q0=q0) #
    dbr = prolif_d(br, d1, q1, ad1, aq1, x0, q0=q0) #n
    before = prolif_d(t, d1, q1, ad1, aq1, x0, q0=q0)
    after = prolif_d(t - br, d2, q2, ad2, aq2, dbr, q0=qbr)
    answer = before
    if t > br:
        answer=after
    return answer
    
def split_fq(t, d1, d2, q1, q2, ad1, aq1, ad2, aq2, x0, q0=0, br = 20) :
    return split_q(t, d1, d2, q1, q2, ad1, aq1, ad2, aq2,x0, q0=q0, br = 20)/split_tot(t, d1, d2, q1, q2, ad1, aq1, ad2, aq2,x0, q0=q0, br = 20)

def split_fd(t, d1, d2, q1, q2, ad1, aq1, ad2, aq2, x0, q0=0, br = 20) :
    return split_d(t, d1, d2, q1, q2, ad1, aq1, ad2, aq2,x0, q0=q0, br = 20)/split_tot(t, d1, d2, q1, q2, ad1, aq1, ad2, aq2,x0, q0=q0, br = 20)
    

def death_rate(Time_h, card_total, card_total_e16000, v) :
    X = np.array([Time_h]).T
    Y= np.array([card_total]).T
    
    reg = LinearRegression().fit(X, Y)
    reg.score(X, Y)
    
    untreated_division_rate = reg.coef_ * np.log(2)
    half_imt = np.log(2)/untreated_division_rate/2
    
    y = np.array(card_total_e16000) - card_total_e16000[0]
    
    def func(a):
        return np.array([np.log2(split_tot(t, untreated_division_rate, v['d'], 0, v['q'], 0, 0, a, a, 1, q0=0, br = half_imt))[0][0] for t in Time_h]) - y
        
    res = least_squares(func, x0=0.003, bounds = (0,np.inf)) 
    
    return half_imt, untreated_division_rate, res["x"]



def death_rate_without_delay(Time_h, card_total, v) :
    X = np.array([Time_h]).T
    Y= np.array([card_total]).T
    
    reg = LinearRegression().fit(X, Y)
    reg.score(X, Y)
    
    untreated_division_rate = reg.coef_ * np.log(2)
    half_imt = np.log(2)/untreated_division_rate/2
    
    y = np.array(card_total) - card_total[0]
    
    def func(a):
        return np.array([np.log2(split_tot(t, v['d'], v['d'], v['q'], v['q'], a, a, a, a, 1, q0=0, br = 0))[0] for t in Time_h]) - y
        
    res = least_squares(func, x0=0.003, bounds = (0,np.inf)) 
    
    return half_imt, untreated_division_rate, res["x"]


######################
#test
"""
mitotic_lifespans = [17.0, 21.7, 21.7, 23.0, 22.5, 20.0, 20.0, 24.2, 19.2, 55.3, 20.7, 29.0, 21.0, 23.7, 34.5, 34.5, 35.3, 35.3, 29.0, 28.3, 57.5, 30.3, 31.2, 23.2, 28.7, 32.8, 21.7, 25.8, 21.3, 61.8, 28.0, 28.7, 22.2, 22.2, 34.0, 27.7, 21.7, 21.7, 28.0, 31.8, 25.0, 20.7, 28.8, 31.7, 45.0, 25.2, 23.2, 42.5, 23.0, 31.8, 18.2, 37.2, 22.2, 22.2, 17.8, 17.8, 41.0, 45.0, 24.2, 60.0]

censored_lifespans = [92.2, 91.5, 91.5, 90.8, 90.8, 90.8, 90.3, 90.3, 90.3, 89.7, 89.7, 89.3, 89.3, 88.2, 87.5, 87.5, 87.3, 86.8, 85.8, 85.8, 84.8, 84.8, 84.2, 84.2, 82.5, 82.3, 82.3, 82.2, 82.2, 82.0, 82.0, 81.5, 81.3, 81.3, 81.5, 80.5, 79.8, 79.8, 79.8, 79.8, 78.7, 78.7, 78.5, 78.0, 77.8, 77.8, 77.8, 77.8, 77.7, 77.7, 77.5, 77.5, 77.5, 77.2, 77.2, 76.7, 76.7, 75.8, 75.8, 75.7, 75.7, 75.2, 74.8, 74.7, 74.3, 74.3, 73.8, 73.8, 73.2, 72.7, 72.7, 72.7, 72.7, 72.2, 72.2, 71.3, 71.3, 71.3, 71.3, 71.3, 71.3, 71.2, 71.2, 71.0, 71.0, 69.7, 69.7]

card_total = [7.73, 7.83, 7.84, 7.94, 7.99, 8.00, 8.09, 8.14, 8.16, 8.22, 8.22, 8.35, 8.36, 8.43, 8.49, 8.55, 8.65, 8.73, 8.80, 8.85, 8.89, 8.94, 8.99, 9.04, 9.09, 9.16, 9.18, 9.29, 9.38, 9.43, 9.51, 9.60, 9.63, 9.69, 9.72, 9.77, 9.80, 9.86, 9.88, 9.92, 9.98, 10.01, 10.07, 10.12, 10.18, 10.23, 10.30, 10.38, 10.43, 10.49, 10.53, 10.57, 10.62, 10.64, 10.68, 10.72, 10.76, 10.81, 10.88, 10.92, 10.96, 11.02, 11.09, 11.14, 11.19, 11.24, 11.30, 11.36, 11.41, 11.47, 11.51, 11.55, 11.59]

card_total_e16000 = [6.86, 6.92, 6.99, 7.00, 7.15, 7.15, 7.23, 7.28, 7.29, 7.27, 7.35, 7.32, 7.38, 7.35, 7.36, 7.48, 7.48, 7.53, 7.53, 7.57, 7.58, 7.54, 7.56, 7.60, 7.69, 7.72, 7.73, 7.75, 7.75, 7.79, 7.79, 7.85, 7.84, 7.83, 7.85, 7.88, 7.89, 7.89, 7.94, 7.91, 7.92, 7.94, 7.95, 7.97, 7.98, 7.95, 7.97, 7.97, 7.97, 7.98, 8.00, 7.97, 8.01, 7.98, 7.97, 8.01, 8.02, 8.04, 8.04, 8.06, 8.07, 8.06, 7.99, 8.03, 8.07, 8.06, 8.03, 8.02, 8.01, 8.04, 8.02, 8.04, 8.05]

Time_h = range(len(card_total))

plt.plot(Time_h, card_total)



Q=0.5

m = mle(mitotic_lifespans, censored_lifespans)
m['x']

#         mu       sigma      lambda           Q 
#18.89452189  1.67126329  0.09735639  0.59049561 


v = q_rates(m)

#         d          q
#0.02475328 0.03569364 

half_imt, untreated_division_rate, a = death_rate(Time_h, card_total, card_total_e16000, v)

#          a 
#0.004111099 


y = np.array(card_total_e16000) - card_total_e16000[0]

 
model_total = np.array([np.log2(split_tot(t, untreated_division_rate, v['d'], 0.0, v['q'], 0, 0, a, a, 1, 0, half_imt))[0][0] for t in Time_h])

model_fq = np.array([split_fq(t, untreated_division_rate, v['d'], 0.0, v['q'], 0, 0, a, a, 1, 0, half_imt)[0][0] for t in Time_h]) * model_total

model_fd = np.array([split_fd(t, untreated_division_rate, v['d'], 0.0, v['q'], 0, 0, a, a, 1, 0, half_imt)[0][0] for t in Time_h]) * model_total



plt.scatter(Time_h, y, s=40, facecolors='none', edgecolors='black')
plt.plot(Time_h, model_total, c="green")
plt.plot(Time_h, model_fq, c="red")
plt.plot(Time_h, model_fd, c="blue")
"""
"""

x0 = np.array([np.median(mitotic_lifespans), np.std(mitotic_lifespans), 1/np.mean(mitotic_lifespans), 0.5])

def qsurvival_nllik(x0):
    mu, sigma, lambd, Q = x0   
    #global mitotic_lifespans, censored_lifespans
    return - sum([np.log(1-Q) + np.log(demg(i, mu = mu, sigma = sigma, lambd = lambd)) for i in  mitotic_lifespans]) - sum([np.log(Q+(1-Q)*pemg(i,  mu = mu, sigma = sigma, lambd = lambd, lower_tail = False)) for i in censored_lifespans])

bounds = [(1,60), (0.1, (max(mitotic_lifespans) - min(mitotic_lifespans))), (0.01, 10), (1e-31, 0.9999999)]
#
#m = minimize(qsurvival_nllik, x0, method='nelder-mead', bounds=bounds) #Change the method to minimize
#x: array([18.89482881,  1.67120341,  0.09735528,  0.59049594])
#x0: array([25.5       , 10.35847627,  0.03450854,  0.5       ])

#m = minimize(qsurvival_nllik, x0, method="L-BFGS-B", bounds=bounds)

#m = minimize(qsurvival_nllik, x0, method="TNC", bounds=bounds)
m = minimize(qsurvival_nllik, x0, method="SLSQP", bounds=bounds)
""" 
    
