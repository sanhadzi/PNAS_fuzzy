#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  2 14:05:54 2019

@author: sanhadzi
"""

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  1 12:42:35 2019

@author: sanhadzi
"""
from lmfit import minimize, Parameters, report_fit
from numpy import loadtxt, savetxt, array, transpose, hsplit, log, exp, diff, delete, matrix, arange, linspace, rint, c_, column_stack
from scipy import stats
from scipy.misc import factorial
import pylab
import matplotlib.pyplot as plt
from sympy import Function, symbols, lambdify, diff, simplify, Matrix, expand
from numpy import arange, exp, array, log, zeros, insert, where, copy
import matplotlib.pyplot as plt
from lmfit import minimize, Parameters, report_fit
from scipy.optimize import brentq, fsolve
import pandas as pd
import seaborn as sns
import numpy as np


# D-, E-, H0, K+, R+
letter_val = {'A': 1.61,'C': 0.35,'D': 0.31,'E': 0.45,'F': 0.29,'G': 0.05,'H': 0.38,'I': 0.44,'K': 0.82,'L': 0.96,'M':0.63,'N': 0.31,'P': 0.01,'Q': 0.56,'R': 1.2,'S': 0.38,'T': 0.14,'V': 0.23,'W': 0.35,'Y': 0.45}        # chakrabatky22
letter_val.update((x, y*2.4) for x, y in letter_val.items())


def w_values(s):
    res = 0
    for char in s:
        if char in letter_val:
            res += letter_val[char]
    return res

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth


sekvenca10='SDLACRLLGQSMDESGLPQLTSYDCEVNAPIQGSRNLLQGEELLRALDQVN' #776-826;
kk=677
#>>>>>>>>>>>>>>>>PPcheck HOTSPOT RESIDUES based on PP check
hotspots=array([(102+kk),(118+kk),(127+kk),(135+kk),(141+kk),(145+kk),(148+kk)]) 

sekvenca=sekvenca10
res, S2 = loadtxt('/???????/hifA_2_nmr_data.txt', unpack=True)#, skiprows=0)


x= symbols("x")
hotspots=hotspots
x= symbols("x")
for i in hotspots:
    print sekvenca[i-775-1],i
    
    
N = len(sekvenca)
N2=len(sekvenca)
print 'IDP length', len(sekvenca)
v=0.048
print hotspots
interakcije=[] ; j=775
for i in sekvenca:
    j=j+1
    print 'ak', i, j
    if any(hotspots==(j)):
        print '>',i,j, 'tuuuu'
        interakcije.append(1)
    else:
        interakcije.append(0)

print 'intractions', interakcije        
  
produkt=1 ; st_odv_x=[]
for i, g in zip(sekvenca,interakcije):
    w=w_values(i)    
    if g==0:
        M= matrix( [[w,v,0], [0,0,1], [v,v,1]])
    if g==1:
        M= Matrix( [[w*x,v,0], [0,0,1], [v,v,1]]) 
        st_odv_x.append(1)
    produkt=M*produkt
    
Q13= Matrix( [[0,0,1]] ) ;     COL_vektor3= Matrix( [[0], [1], [1]])
Z3= Q13*produkt*COL_vektor3 ; Z3=Z3[0]
NOx=sum(st_odv_x) 

dZ3=diff(Z3,x, NOx)
Z3bound=((dZ3*x**(NOx))/factorial(NOx))


print 'print partition function', expand(Z3bound)
 
lam_Z3bound=lambdify((x), Z3bound, modules='numpy')

j=0
probab=[] ; sekv=[]
produkt=1 ; produkt1=1 ; produkt0=1 ; xx=0.5

for i,g in zip(sekvenca,interakcije):
    j=j+1
    w_i=w_values(i)
    if g==0:        
        M3w= matrix( [[w_i,0,0], [0,0,0], [0,0,0]])
    if g==1:        
        M3w= matrix( [[w_i*x,0,0], [0,0,0], [0,0,0]])   
        
    produkt=produkt*M3w
    #pohodni ostanki
    for i,g in (zip(sekvenca,interakcije))[j:]:
        w1=w_values(i)   
        if g==0:        
            M1= matrix( [[w1,v,0], [0,0,1], [v,v,1]])   
        if g==1:        
            M1= matrix( [[w1*x,v,0], [0,0,1], [v,v,1]])    
        produkt1=produkt1*M1

    #predhodni ostanki        
    for i,g in (zip(sekvenca,interakcije))[:(j-1)]:
        w0=w_values(i)
        if g==0:        
            M0= matrix( [[w0,v,0], [0,0,1], [v,v,1]])   
        if g==1:        
            M0= matrix( [[w0*x,v,0], [0,0,1], [v,v,1]]) 
        produkt0=produkt0*M0
    
    pp=produkt0*produkt*produkt1        
    Z_i=Q13*pp*COL_vektor3 
    dZ_i=diff(Z_i,x, NOx)
    Z3_i_B=((dZ_i*x**(NOx))/factorial(NOx))
    Z3_i_B=Z3_i_B[0]
    lam_Z3_i_B=lambdify((x), Z3_i_B, modules='numpy')
    
    ZZ3_i_B=lam_Z3_i_B(xx) ; ZZbound=lam_Z3bound(xx)

    p2=ZZ3_i_B/ZZbound #; p2=p2.flat[0]
    produkt=1 ; produkt1=1 ; produkt0=1
    probab.append(p2)
    sekv.append(j+776)        
    

prob2=smooth(probab,3)
    
print probab
print prob2
print sekv

for i,j,k,l in zip(sekv,res,prob2,S2):
    print i,j,'__',k,l
    
#fig = plt.figure()
#plt.plot(sekv, probab, 'g-')
plt.plot(sekv, prob2, 'k-')
plt.bar(res, S2, )
