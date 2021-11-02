#!/usr/bin/python
# -*- coding: utf-8 -*-
import numpy as np
from numpy import linalg as LA
from numpy.linalg import matrix_power
import time 
import os 
import datetime 
import sys
import random
import math
import scipy as sp
import scipy.linalg
from scipy.linalg import expm
from matplotlib.pyplot import *
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot

startTime = time.time()
print ("STARTED:" , datetime.datetime.now().strftime("%d %B %Y %H:%M:%S"))

if len(sys.argv) < 7:
  print("Usage: python", str(sys.argv[0]), "READ-IN? " "SAVE-or-NOT? " "NCOL " "NITERS " "D " "LAMBDA ")
  sys.exit(1)

READIN = int(sys.argv[1])
SAVE = int(sys.argv[2])
NCOL = int(sys.argv[3]) 
Niters_sim = int(sys.argv[4]) 
NMAT = int(sys.argv[5])
LAMBDA = float(sys.argv[6])
if NMAT < 2:
    print ("NMAT must be at least two")
    sys.exit(1) 
if READIN not in [0,1]:
    print ("Wrong input for READIN")
    sys.exit(1)
if SAVE not in [0,1]:
    print ("Wrong input for SAVE")
    sys.exit(1)

COUPLING = float(NCOL/(4.0*LAMBDA))
GENS = NCOL**2 - 1
dt = 5e-4
nsteps = int(1e-2/dt)
GAP = 1
t2 = np.zeros((NMAT),dtype=float)
t4 = np.zeros((NMAT),dtype=float)
X = np.zeros((NMAT, NCOL, NCOL), dtype=complex)
mom_X = np.zeros((NMAT, NCOL, NCOL), dtype=complex)
f_X = np.zeros((NMAT, NCOL, NCOL), dtype=complex)
X_bak = np.zeros((NMAT, NCOL, NCOL), dtype=complex)
HAM, expDH, ACT, scalar = [],[],[],[]

print ("Yang-Mills type matrix model with %2.0f matrices" % (NMAT)) 
print ("NCOL = " "%3.0f " ","  " and coupling = " " %4.2f" % (NCOL, COUPLING)) 
print ("--------------------------------------------")

def dagger(a):
    return np.transpose(a).conj()

def box_muller():
    PI = 2.0*math.asin(1.0);    
    r = random.uniform(0,1)
    s = random.uniform(0,1)
    p = np.sqrt(-2.0*np.log(r)) * math.sin(2.0*PI*s)
    q = np.sqrt(-2.0*np.log(r)) * math.cos(2.0*PI*s)
    return p,q

def comm(A,B):
    return np.dot(A,B) - np.dot(B,A)

def unit_matrix():
    matrix = np.zeros((NCOL, NCOL), dtype=complex)
    for i in range (NCOL):
        matrix[i][i] = complex(1.0,0.0)
    return matrix

def copy_fields(b):
    for j in range(NMAT):
        X_bak[j] = b[j]
    return X_bak

def rejected_go_back_old_fields(a):
    for j in range(NMAT):
        X[j] = a[j]
    return X

def refresh_mom():
    for j in range (NMAT):
        mom_X[j] = random_hermitian()
    return mom_X

def random_hermitian():
    tmp = np.zeros((NCOL, NCOL), dtype=complex)
    for i in range (NCOL):

        for j in range (i+1, NCOL):
            r1, r2 = box_muller()
            tmp[i][j] = complex(r1, r2)/math.sqrt(2)
            tmp[j][i] = complex(r1, -r2)/math.sqrt(2)

    for i in range (NCOL):
        r1, r2 = box_muller()
        tmp[i][i] = complex(r1, 0.0)
    return tmp 

def makeH(tmp):
    tmp2 = 0.50*(tmp+dagger(tmp)) - (0.50*np.trace(tmp+dagger(tmp))*np.eye(NCOL))/NCOL
    for i in range (NCOL):
        tmp2[i][i] = complex(tmp[i][i].real,0.0)  
    if np.allclose(tmp2, dagger(tmp2)) == False:
        print ("WARNING: Couldn't make hermitian.")
    return tmp2

def hamil(mom_X):
    s = 0.0 
    for j in range (NMAT):
        s += 0.50 * np.trace(np.dot(dagger(mom_X[j]),mom_X[j])).real
    return s    

def potential(X):
    s1 = 0.0 
    for i in range (NMAT):
        for j in range (i+1, NMAT): 
            co = np.dot(X[i],X[j]) - np.dot(X[j],X[i])
            tr = np.trace(np.dot(co,co))
            s1 -= COUPLING*tr.real 
    return s1

def force(X):

    tmp_X = np.zeros((NMAT, NCOL, NCOL), dtype=complex)
    for i in range (NMAT): 
        for j in range (NMAT):
            if i == j:
                continue 
            else:
                temp = comm(X[i], X[j])
                tmp_X[i] -= comm(X[j], temp)
        f_X[i] = 2.0*COUPLING*dagger(tmp_X[i])

    for j in range(NMAT):
        if np.allclose(f_X[j], dagger(f_X[j])) == False:
            f_X[j] = makeH(f_X[j])
    return f_X  

def leapfrog(X,mom_X, dt):
    for j in range(NMAT):
        X[j] += mom_X[j] * dt/2.0
    f_X = force(X)

    for step in range(nsteps):
        for j in range(NMAT):
            mom_X[j] -= f_X[j] * dt
            X[j] += mom_X[j] * dt
        f_X = force(X)

    for j in range(NMAT):
        mom_X[j] -= f_X[j] * dt
        X[j] += mom_X[j] * dt/2.0
    
    return X, mom_X, f_X

def update(X):
    mom_X = refresh_mom()
    s1 = hamil(mom_X)
    s2 = potential(X)
    start_act =  s1 + s2
    X_bak = copy_fields(X) 
    X, mom_X, f_X = leapfrog(X,mom_X,dt)
    s1 = hamil(mom_X)
    s2 = potential(X)
    end_act = s1 + s2
    change = end_act - start_act
    HAM.append(abs(change))
    expDH.append(np.exp(-1.0*change))   

    if np.exp(-change) < random.uniform(0,1):
        X = rejected_go_back_old_fields(X_bak)
        print(("REJECT: deltaH = " "%10.7f " " startH = " "%10.7f" " endH = " "%10.7f" % (change, start_act, end_act)))
    else:
        print(("ACCEPT: deltaH = " "%10.7f " "startH = " "%10.7f" " endH = " "%10.7f" % (change, start_act, end_act)))

    ACT.append(s2)
    tmp = 0.0 
    for i in range (0,NMAT):
        val = np.trace(X[i] @ X[i]).real/NCOL
        val2 = np.trace(X[i] @ X[i] @ X[i] @ X[i]).real/NCOL
        t2[i] = val 
        t4[i] = val2 
        tmp += val 

    tmp /= NMAT 
    scalar.append(tmp) 
    if MDTU%GAP == 0:
        f3.write("%4.8f \n" % (s2/GENS))
        for item in t2:
            f4.write("%4.8f " % item)
        for item in t4:
            f5.write("%4.8f " % item)
        f4.write("\n")
        f5.write("\n") 

    return X


if __name__ == '__main__':

    if READIN == 0:
        print ("Starting from fresh")
        for i in range (NMAT):  
            X[i] = 0.0  

    if READIN == 1:
        name_f = "config_YM_N{}_l_{}_D_{}.npy".format(NCOL, LAMBDA, NMAT)

        if os.path.isfile(name_f) == True: 
            print ("Reading old configuration file:", name_f)
            
            A = np.load(name_f)
            for i in range (NMAT):
                for j in range (NCOL):
                    for k in range (NCOL):
                        X[i][j][k] = A[i][j][k] 

            for j in range(NMAT):
                if np.allclose(X[j], dagger(X[j])) == False:
                    print ("Input configuration not hermitian, making it so")
                    X[j] = makeH(X[j])
        else: 
            print ("Can't find config. file for this NCOL and LAM")
            print ("Starting from fresh")
            for i in range (NMAT):  
                X[i] = 0.0

    f3 = open('action_N%s_D%s.txt' %(NCOL,NMAT), 'w')
    f4 = open('t2_N%s_D%s.txt' %(NCOL,NMAT), 'w')
    f5 = open('t4_N%s_D%s.txt' %(NCOL,NMAT), 'w')

    for MDTU in range (1, Niters_sim+1): 
        X = update(X)

        if MDTU%10 == 0 and SAVE == 1:
            name_f = "config_YM_N{}_l_{}_D_{}.npy".format(NCOL, LAMBDA, NMAT)
            print ("Saving configuration file: ", name_f)
            np.save(name_f, X)

    ACT = [x/GENS for x in ACT] 
    f3.close()
    f4.close()
    f5.close()
    print ("--------------------------------------------")
    print("<S> = ", np.mean(ACT), "+/-", (np.std(ACT)/np.sqrt(np.size(ACT) - 1.0)))
    print ("COMPLETED:" , datetime.datetime.now().strftime("%d %B %Y %H:%M:%S")) 
    endTime = time.time() 
    # Plot results!
    t2plot = plt.figure(1) 
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    MDTU = np.linspace(0, int(Niters_sim/GAP), int(Niters_sim/GAP), endpoint=True)
    plt.ylabel(r'$\langle R^2 \rangle$')
    plt.xlabel('Time units')
    plot(MDTU, scalar, 'teal') 
    plt.grid(which='major', axis='y', linestyle='--')
    act_plot = plt.figure(2) 
    plt.ylabel(r'$\langle S/(N^2-1) \rangle$')
    plt.xlabel('Time units')
    plt.axhline(y = NMAT/4.0, color='blue', linestyle='--')
    plot(MDTU, ACT, 'blue') 
    plt.grid(which='major', axis='y', linestyle='--')
    outname = "YM_N%s_D%s" %(NCOL,NMAT)
    pp = PdfPages(outname+'.pdf')
    pp.savefig(t2plot, dpi = 300, transparent = True)
    pp.savefig(act_plot, dpi = 300, transparent = True)
    pp.close()
    print ("Running time:", round(endTime - startTime, 2),  "seconds")