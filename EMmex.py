import msprime
import tskit
import random
import numpy as np
import random
import sklearn

import pandas as pd
from numpy import linalg as LNG 
from random import randint, randrange


from sklearn.metrics import confusion_matrix
from sklearn.metrics import classification_report 


import seaborn
import matplotlib
from matplotlib import pyplot as plt

import math
import scipy
import scipy.optimize as optimize
from scipy.optimize import minimize



import HMMmex as hmm

def point_in_set(p, m_a):
    if m_a==[[]]:
        return False
    def point_in(p, a):
        if p>=a[0] and p<=a[1]:
            return True
        else:
            return False
    f=False
    for j in m_a:
        f=point_in(p,j)
        if f==True:
            return f
    return f



def intersections(a,b):
    ranges = []
    i = j = 0
    while i < len(a) and j < len(b):
        a_left, a_right = a[i]
        b_left, b_right = b[j]

        if a_right < b_right:
            i += 1
        else:
            j += 1

        if a_right >= b_left and b_right >= a_left:
            end_pts = sorted([a_left, a_right, b_left, b_right])
            middle = [end_pts[1], end_pts[2]]
            ranges.append(middle)

    ri = 0
    while ri < len(ranges)-1:
        if ranges[ri][1] == ranges[ri+1][0]:
            ranges[ri:ri+2] = [[ranges[ri][0], ranges[ri+1][1]]]
        ri += 1
    return ranges

N=5

# forward-algo
def alpha_scaled_opt(a,b, o, p):
    
    c = np.zeros(len(o)) #scaling factors, которые как раз позволяют не обнулиться
    
    alpha = np.zeros((N, len(o)))
    alpha[:, 0] = b[:, o[0][0],o[0][1],o[0][2],o[0][3]] * p
    
    c[0] = 1 / sum(alpha[:, 0])
    alpha[:, 0] = alpha[:, 0] / sum(alpha[:, 0])
    for t in range(1, len(o)):   
        
        for i in range(0, N):
            alpha[i, t] = np.dot(alpha[:, t-1],a[:,i]) * b[i, o[t][0],o[t][1],o[t][2],o[t][3]] 
            
        c[t] = 1 / sum(alpha[:,t]) #сохраняем множители        
        alpha[:, t] = alpha[:, t] / sum(alpha[:,t])     
        
    return alpha, c

# Backward procedure. Scaled case.
def beta_scaled_opt(a,b, o, scaling_factors):
    
    beta = np.zeros((N, len(o)))
    
    length = len(o)
    beta[:, length - 1] = np.ones(N)*scaling_factors[length-1] 
    

    for t in range(len(o)-2,-1,-1):             
        for i in range(0, N):             
            for l in range(0, N):
                beta[i, t] += a[i, l] * b[l, o[t+1][0],o[t+1][1],o[t+1][2],o[t+1][3]] * beta[l, t+1]
                
        beta[:, t] = beta[:, t] * scaling_factors[t]

    return beta 

# gamma matrix
def def_gamma(alpha, beta):
        
    gamma = np.zeros((N,len(alpha[0])))
    for m in range(0, len(alpha[0])):
        denom = sum(alpha[:, m]*beta[:,m])
        
        for i in range(0,N):
            gamma[i, m] = (alpha[i, m] * beta[i, m]) / denom
    return gamma


# ksi[i, j, t]
def def_ksi( a, b, o, alpha, beta):
    
    M = len(o)
    ksi = np.zeros((N, N, M-1))
    
    for t in range(0, M-1):
        
        denom = 0
        for i in range(0, N):
            for j in range(0, N):
                denom += alpha[i, t] * a[i, j] * b[j, o[t+1][0],o[t+1][1],o[t+1][2],o[t+1][3]] * beta[j, t+1]
                
        
        for i in range(0, N):
            for j in range(0, N):
                ksi[i, j, t] = (alpha[i, t]*a[i, j]*b[j, o[t+1][0],o[t+1][1],o[t+1][2],o[t+1][3]] * beta[j, t+1]) / denom
    
    return ksi




def new_lambda_mex(o, gamma):
    lmbd=0
    
    for t in range(1, len(o), 1):
        lmbd += o[t, 0] * ( gamma[0, t] + gamma[1, t]) + o[t, 1] * (gamma[2, t] + gamma[3, t]) + o[t, 2]  * gamma[4, t]     
        
    return lmbd/(len(o)-1)


def new_lambda_n(o, gamma):
    lmbd = 0
    
    for t in range(1, len(o), 1):
        lmbd += o[t, 3] * ( gamma[0, t] + gamma[2, t]+gamma[4, t]) + o[t, 2] * (gamma[1, t] + gamma[3, t]) 
        
    return lmbd/ (len(o)-1)

def new_lambda_i(o, gamma):
    nom=0
    denom=0
    
    for t in range(1, len(o), 1):   
        nom += o[t, 3] * (gamma[1, t] + gamma[3, t])
        denom += gamma[1, t] + gamma[3, t]
        
    return nom / denom

def new_lambda_af(o, gamma):
    nom, denom = 0, 0
    
    for t in range(1, len(o), 1):   
        nom += o[t, 2] * ( gamma[0, t] + gamma[2, t]) + (o[t, 0] + o[t, 1]) * gamma[4, t]
        denom += gamma[0, t] + gamma[2, t]+ 2 * gamma[4, t]    
    
    return nom/ denom

def new_lambda_ea(o, gamma):
    nom, denom = 0, 0
    for t in range(1, len(o), 1):  
        nom += o[t, 1] * ( gamma[0, t] + gamma[1, t]) + o[t, 0]  * (gamma[2, t]+gamma[3,t])
        denom += gamma[0, t] + gamma[1, t] +  gamma[2, t] + gamma[3, t]
    
    return nom/ denom



def new_a(coeff_a, lmbd_0, bnds, cut, mut_rate,rr):
    d=cut* mut_rate
    x0=[lmbd_0[5], lmbd_0[6]]
    def multi_Q(x):
        x=np.array(x)
        Q=0            
        a = hmm.initA(x[0]/d, x[1]/d, rr, cut, lmbd_0[7], lmbd_0[8], lmbd_0[9], lmbd_0[10])            



        for ii in range(N):
            for jj in range(N):
                Q += math.log(a[ii,jj]) * coeff_a[ii,jj]

        return -Q

    def gradient_respecting_bounds(bounds, fun, eps=1e-8):
        """bounds: list of tuples (lower, upper)"""
        def gradient(x):
            fx = fun(x)
            grad = np.zeros(len(x))
            for k in range(len(x)):
                d = np.zeros(len(x))
                d[k] = eps if x[k] + eps <= bounds[k][1] else -eps
                grad[k] = (fun(x + d) - fx) / d[k]
            return grad
        return gradient

    opt_result = scipy.optimize.minimize(multi_Q, x0,  bounds=bnds,  
                                         jac=gradient_respecting_bounds(bnds, multi_Q),  
                                         method="L-BFGS-B" )


    return opt_result.x[0], opt_result.x[1]            


        
        
    

def E_step(cut, p, o, n_states, mut_rate, rr,lambda_old, bnds):

    d=mut_rate * cut
    b = hmm.initB(mut_rate, cut, lambda_old[0:5], n_states)  
    a = hmm.initA(lambda_old[5]/d, lambda_old[6]/d, rr, cut, lambda_old[7], lambda_old[8], lambda_old[9], lambda_old[10])
    
    alpha, sc_factors = alpha_scaled_opt(a,b, o, p)
    beta = beta_scaled_opt(a, b, o, sc_factors)    
    gamma = def_gamma(alpha, beta)
    ks = def_ksi( a, b, o, alpha, beta)
    
    
    coeff_a = np.zeros((N, N))
    for i in range(N):
        for j in range(N):                
            for t in range(len(o)-1):
                coeff_a[i, j] += ks[i, j, t] 
    
    upd = new_a(coeff_a, lambda_old, bnds, cut, mut_rate,rr)

    
    return [new_lambda_i(o, gamma), new_lambda_n(o, gamma), new_lambda_af(o, gamma), 
            new_lambda_ea(o, gamma),new_lambda_mex(o, gamma), upd[0], upd[1], 
            lambda_old[7], lambda_old[8], lambda_old[9], lambda_old[10]]


def EM_algorithm(p, o, n_states, mut_rate, rr, lambda_0, epsilon, cut, bnds):     

    lmbd = np.array(lambda_0)
    em_steps = 20

    for i in range(em_steps):
        lmbd_new = np.array(E_step(cut, p, o, n_states, mut_rate, rr,lambda_0, bnds))


        if LNG.norm(lmbd_new-lmbd) < epsilon:
            break
        lmbd = lmbd_new

    return lmbd







#EM where we use several observations simultaneously
def new_lambda_mex_common(o_mas, gamma_mas):
    lmbd=0
    for i in range(len(o_mas)):
        o=o_mas[i]
        gamma=gamma_mas[i]
        for t in range(1, len(o), 1):
            lmbd += o[t, 0] * ( gamma[0, t] + gamma[1, t]) + o[t, 1] * (gamma[2, t] + gamma[3, t]) + o[t, 2]  * gamma[4, t]     
        
    return lmbd/((len(o)-1)*len(o_mas))


def new_lambda_n_common(o_mas, gamma_mas):
    lmbd = 0
    for i in range(len(o_mas)):
        o=o_mas[i]
        gamma=gamma_mas[i]
        for t in range(1, len(o), 1):
            lmbd += o[t, 3] * ( gamma[0, t] + gamma[2, t]+gamma[4, t]) + o[t, 2] * (gamma[1, t] + gamma[3, t]) 
        
    return lmbd/ ((len(o)-1)*len(o_mas))

def new_lambda_i_common(o_mas, gamma_mas):
    nom=0
    denom=0
    for i in range(len(o_mas)):
        o=o_mas[i]
        gamma=gamma_mas[i]   
        for t in range(1, len(o), 1):   
            nom += o[t, 3] * (gamma[1, t] + gamma[3, t])
            denom += gamma[1, t] + gamma[3, t]
        
    return nom / denom

def new_lambda_af_common(o_mas, gamma_mas):
    nom, denom = 0, 0
    for i in range(len(o_mas)):
        o=o_mas[i]
        gamma=gamma_mas[i] 
        for t in range(1, len(o), 1):   
            nom += o[t, 2] * ( gamma[0, t] + gamma[2, t]) + (o[t, 0] + o[t, 1]) * gamma[4, t]
            denom += gamma[0, t] + gamma[2, t]+ 2 * gamma[4, t]    
    
    return nom/ denom

def new_lambda_ea_common(o_mas, gamma_mas):
    nom, denom = 0, 0
    for i in range(len(o_mas)):
        o=o_mas[i]
        gamma=gamma_mas[i]
        for t in range(1, len(o), 1):  
            nom += o[t, 1] * ( gamma[0, t] + gamma[1, t]) + o[t, 0]  * (gamma[2, t]+gamma[3,t])
            denom += gamma[0, t] + gamma[1, t] +  gamma[2, t] + gamma[3, t]
    
    return nom/ denom

def new_a_common(coeff_a, lmbd_0, bnds, cut, mut_rate,rr):
    d=cut* mut_rate
    x0=[lmbd_0[5], lmbd_0[6]]
    def multi_Q(x):
        x=np.array(x)
        Q=0            
        a = hmm.initA(x[0]/d, x[1]/d, rr, cut, lmbd_0[7], lmbd_0[8], lmbd_0[9], lmbd_0[10])            



        for ii in range(N):
            for jj in range(N):
                Q += math.log(a[ii,jj]) * coeff_a[ii,jj]

        return -Q

    def gradient_respecting_bounds(bounds, fun, eps=1e-8):
        """bounds: list of tuples (lower, upper)"""
        def gradient(x):
            fx = fun(x)
            grad = np.zeros(len(x))
            for k in range(len(x)):
                d = np.zeros(len(x))
                d[k] = eps if x[k] + eps <= bounds[k][1] else -eps
                grad[k] = (fun(x + d) - fx) / d[k]
            return grad
        return gradient

    opt_result = scipy.optimize.minimize(multi_Q, x0,  bounds=bnds,  
                                         jac=gradient_respecting_bounds(bnds, multi_Q),  
                                         method="L-BFGS-B" )


    return opt_result.x[0], opt_result.x[1]       



def E_step_common(cut, p, o_mas, n_states, mut_rate, rr,lambda_old, bnds):

    d=mut_rate * cut
    gamma_mas, ks_mas=[], []    
    coeff_a = np.zeros((N, N))

    b = hmm.initB(mut_rate, cut, lambda_old[0:5], n_states)  
    a = hmm.initA(lambda_old[5]/d, lambda_old[6]/d, rr, cut, lambda_old[7], lambda_old[8], lambda_old[9], lambda_old[10])
    for o in o_mas:
        
        alpha, sc_factors = alpha_scaled_opt(a,b, o, p)
        beta = beta_scaled_opt(a, b, o, sc_factors)    
        gamma = def_gamma(alpha, beta)
        ks = def_ksi( a, b, o, alpha, beta)

        gamma_mas.append(gamma)
        ks_mas.append(ks)      
        
    

        for i in range(N):
            for j in range(N):                
                for t in range(len(o)-1):
                    coeff_a[i, j] += ks[i, j, t] 
    
    upd = new_a(coeff_a, lambda_old, bnds, cut, mut_rate,rr)

    
    return [new_lambda_i_common(o_mas, gamma_mas), new_lambda_n_common(o_mas, gamma_mas), 
            new_lambda_af_common(o_mas, gamma_mas), 
            new_lambda_ea_common(o_mas, gamma_mas),new_lambda_mex_common(o_mas, gamma_mas), 
            upd[0], upd[1], 
            lambda_old[7], lambda_old[8], lambda_old[9], lambda_old[10]]



def EM_common(p, o_mas, n_states, mut_rate, rr, lambda_0, epsilon, cut, bnds, em_steps):
    d=mut_rate * cut
    lmbd = np.array(lambda_0)


    for i in range(em_steps):
#        print('Step ', i, lmbd[0:7]/d*29)
        lmbd_new = np.array(E_step_common(cut, p, o_mas, n_states, mut_rate, rr,lmbd, bnds))


        if LNG.norm(lmbd_new-lmbd) < epsilon:
            break
        lmbd = lmbd_new

    return lmbd
 
 
 
 
 
 
 
 
 
 
 
 
 
    
a_gaps=np.identity(5)
b_gaps=np.ones((5,1))

    
    
#EM-common+gaps
 # forward-algo
def alpha_scaled_opt_gaps(a,b, o, p, gaps):   
    c = np.zeros(len(o)) #scaling factors, которые как раз позволяют не обнулиться
    
    alpha = np.zeros((N, len(o)))
    alpha[:, 0] = b[:, o[0][0],o[0][1],o[0][2],o[0][3]] * p
    
    c[0] = 1 / sum(alpha[:, 0])
    alpha[:, 0] = alpha[:, 0] / sum(alpha[:, 0])
    for t in range(1, len(o)):   
        if point_in_set(t, gaps)==True:

            for i in range(0, N):
                
                alpha[i, t] = np.dot(alpha[:, t-1],a_gaps[:,i]) * b_gaps[i][0]     

        else:
        
            for i in range(0, N):
            
                alpha[i, t] = np.dot(alpha[:, t-1],a[:,i]) * b[i, o[t][0],o[t][1],o[t][2],o[t][3]] 
            
        c[t] = 1 / sum(alpha[:,t]) #сохраняем множители        
        alpha[:, t] = alpha[:, t] / sum(alpha[:,t])     
        
    return alpha, c

# Backward procedure. Scaled case.
def beta_scaled_opt_gaps(a,b, o, scaling_factors, gaps):
    
    beta = np.zeros((N, len(o)))
    
    length = len(o)
    beta[:, length - 1] = np.ones(N)*scaling_factors[length-1] 
    

    for t in range(len(o)-2,-1,-1):        
        if point_in_set(t, gaps)==True:
            for i in range(0, N):             
                for l in range(0, N):
                    beta[i, t] += a_gaps[i, l] * b_gaps[l][0] * beta[l, t+1]       
        else:     
            for i in range(0, N):             
                for l in range(0, N):
                    beta[i, t] += a[i, l] * b[l, o[t+1][0],o[t+1][1],o[t+1][2],o[t+1][3]] * beta[l, t+1]
                
        beta[:, t] = beta[:, t] * scaling_factors[t]

    return beta 

# gamma matrix
def def_gamma_gaps(alpha, beta, gaps):
        
    gamma = np.zeros((N,len(alpha[0])))
    for m in range(0, len(alpha[0])):
        denom = sum(alpha[:, m]*beta[:,m])
        
        for i in range(0,N):
            gamma[i, m] = (alpha[i, m] * beta[i, m]) / denom
    return gamma


# ksi[i, j, t]
def def_ksi_gaps( a, b, o, alpha, beta, gaps):
    
    M = len(o)
    ksi = np.zeros((N, N, M-1))
    
    for t in range(0, M-1):
        if point_in_set(t, gaps)==True:
            denom = 0
            for i in range(0, N):
                for j in range(0, N):
                    denom += alpha[i, t] * a_gaps[i, j] * b_gaps[j][0] * beta[j, t+1]        
            for i in range(0, N):
                for j in range(0, N):       
                    ksi[i, j, t] = (alpha[i, t]*a_gaps[i, j]*b_gaps[j,0]* beta[j, t+1]) / denom
        else:
        
            denom = 0
            for i in range(0, N):
                for j in range(0, N):
                    denom += alpha[i, t] * a[i, j] * b[j, o[t+1][0],o[t+1][1],o[t+1][2],o[t+1][3]] * beta[j, t+1]
                
        
            for i in range(0, N):
                for j in range(0, N):
                    ksi[i, j, t] = (alpha[i, t]*a[i, j]*b[j, o[t+1][0],o[t+1][1],o[t+1][2],o[t+1][3]] * beta[j, t+1]) / denom
    
    return ksi

def new_lambda_mex_gaps(o, gamma, gaps):
    lmbd=0
    
    s=0    
    for j in gaps:
        s+=j[1]-j[0]+1
    
    for t in range(1, len(o), 1):
        if point_in_set(t, gaps)==False:
            lmbd += o[t, 0] * ( gamma[0, t] + gamma[1, t]) + o[t, 1] * (gamma[2, t] + gamma[3, t]) + o[t, 2]  * gamma[4, t]     
        
    return lmbd/(len(o)-s-1)


def new_lambda_n(o, gamma,gaps):
    lmbd = 0
    s=0    
    for j in gaps:
        s+=j[1]-j[0]+1
    
    for t in range(1, len(o), 1):
        if point_in_set(t, gaps)==False:
            lmbd += o[t, 3] * ( gamma[0, t] + gamma[2, t]+gamma[4, t]) + o[t, 2] * (gamma[1, t] + gamma[3, t]) 
        
    return lmbd/ (len(o)-s-1)

def new_lambda_i(o, gamma, gaps):
    nom=0
    denom=0
    
    for t in range(1, len(o), 1, gaps):  
        if point_in_set(t, gaps)==False: 
            nom += o[t, 3] * (gamma[1, t] + gamma[3, t])
            denom += gamma[1, t] + gamma[3, t]
        
    return nom / denom

def new_lambda_af(o, gamma, gaps):
    nom, denom = 0, 0
    
    for t in range(1, len(o), 1):   
        if point_in_set(t, gaps)==False:
            nom += o[t, 2] * ( gamma[0, t] + gamma[2, t]) + (o[t, 0] + o[t, 1]) * gamma[4, t]
            denom += gamma[0, t] + gamma[2, t]+ 2 * gamma[4, t]    
    
    return nom/ denom

def new_lambda_ea(o, gamma, gaps):
    nom, denom = 0, 0
    for t in range(1, len(o), 1):  
        if point_in_set(t, gaps)==False:
            nom += o[t, 1] * ( gamma[0, t] + gamma[1, t]) + o[t, 0]  * (gamma[2, t]+gamma[3,t])
            denom += gamma[0, t] + gamma[1, t] +  gamma[2, t] + gamma[3, t]
    
    return nom/ denom 
    
    
    
def E_step_common_gaps(cut, p, o_mas, n_states, mut_rate, rr,lambda_old, bnds, gaps):

    d=mut_rate * cut
    gamma_mas, ks_mas=[], []    
    coeff_a = np.zeros((N, N))

    b = hmm.initB(mut_rate, cut, lambda_old[0:5], n_states)  
    a = hmm.initA(lambda_old[5]/d, lambda_old[6]/d, rr, cut, lambda_old[7], lambda_old[8], lambda_old[9], lambda_old[10])
    for o in o_mas:
        
        alpha, sc_factors = alpha_scaled_opt_gaps(a,b, o, p, gaps)
        beta = beta_scaled_opt_gaps(a, b, o, sc_factors, gaps)    
        gamma = def_gamma_gaps(alpha, beta, gaps)
        ks = def_ksi_gaps( a, b, o, alpha, beta, gaps)

        gamma_mas.append(gamma)
        ks_mas.append(ks)      
        
    

        for i in range(N):
            for j in range(N):                
                for t in range(len(o)-1):
                    if point_in_set(t, gaps)==False:
                        coeff_a[i, j] += ks[i, j, t] 
    
#    upd = new_a(coeff_a, lambda_old, bnds, cut, mut_rate,rr)

    
    return [new_lambda_i_common(o_mas, gamma_mas), new_lambda_n_common(o_mas, gamma_mas), 
            new_lambda_af_common(o_mas, gamma_mas), 
            new_lambda_ea_common(o_mas, gamma_mas),new_lambda_mex_common(o_mas, gamma_mas), 
#            upd[0], upd[1], 
            lambda_old[5], lambda_old[6],
            lambda_old[7], lambda_old[8], lambda_old[9], lambda_old[10]]



def EM_common_gaps(p, o_mas, n_states, mut_rate, rr, lambda_0, epsilon, cut, bnds, em_steps, gaps):
    d=mut_rate * cut
    lmbd = np.array(lambda_0)


    for i in range(em_steps):


        lmbd_new = np.array(E_step_common_gaps(cut, p, o_mas, n_states, mut_rate, rr,lmbd, bnds, gaps))

        print('Step',  LNG.norm(lmbd_new-lmbd))
        if LNG.norm(lmbd_new-lmbd) < epsilon:
            break
        lmbd = lmbd_new

    return lmbd
       
