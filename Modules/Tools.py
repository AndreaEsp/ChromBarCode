# -*- coding: utf-8 -*-
import numpy as np
import scipy
import os,sys

def build_random_model(binding_domain):
    """Randomize a binding domain by bootstrapping"""
    tot = np.sum(binding_domain)
    random_model = np.zeros(binding_domain.shape)
    while tot != 0:
        random_model[np.random.randint(len(random_model))] += 1
        tot -= 1    
    return random_model

def compute_overlap(f_1,f_2):
    """Compute the overlap between a pair of binding domains"""
    if len(f_1)!=len(f_2):
        print('compute_overlap:error dimensions')
    else:
        overlap = np.dot(f_1,f_2)/(np.linalg.norm(f_1)*np.linalg.norm(f_2))
        return overlap

def compute_masked_overlap(f_1,f_2):
    """Same as compute_overlap, but handling masked array"""
    if len(f_1)!=len(f_2):
        print('compute_overlap:error dimensions')
    else:
        overlap = np.ma.dot(f_1,f_2)/( np.ma.sqrt(np.ma.dot(f_1,f_1))*np.ma.sqrt(np.ma.dot(f_2,f_2)) )
        return overlap

def match_by_overlap(polymer0,polymer1):
    """Given two SBS polymers objects, matches their binding domains according to their overlap"""
    max_overlap = np.zeros((polymer0.n,3))
    
    counter = 0
    tot = polymer0.n
    var_0 = np.zeros(polymer0.n)
    var_1 = np.zeros(polymer1.n)
    
    while counter < tot:
        
        over = np.zeros((polymer0.n,polymer1.n))
        for i in np.arange(0,polymer0.n):
            for j in np.arange(0,polymer1.n):
                if(var_0[i]==0 and var_1[j]==0):
                    over[i,j] = compute_overlap(polymer0.body[:,i+1],polymer1.body[:,j+1])
    
        dom_0 = np.where(over==over.max())[0][0]
        dom_1 = np.where(over==over.max())[1][0]
        max_overlap[counter] = [dom_0+1,dom_1+1,over.max()]
        var_0[dom_0]=1
        var_1[dom_1]=1
        counter+=1
    
    return max_overlap


def all_overlaps(polymer0,polymer1):
    """Compute overlap among all the pairs of binding domains of two SBS polymers objects"""
    over = np.zeros((polymer0.n,polymer1.n))
    for i in np.arange(0,polymer0.n):
        for j in np.arange(0,polymer1.n):
            over[i,j] = compute_overlap(polymer0.body[:,i+1],polymer1.body[:,j+1])

    return over


def intra_overlaps(polymer):
    """Compute overlap among all the pairs of binding domains of the given SBS polymer object"""
    over = []
    for i in np.arange(0,polymer.n):
        for j in np.arange(i+1,polymer.n):
            over.append( compute_overlap(polymer.body[:,i+1],polymer.body[:,j+1]) )

    return np.array(over)
