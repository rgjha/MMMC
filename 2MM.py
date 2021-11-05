#!/usr/bin/python
# -*- coding: utf-8 -*- 
import numpy as np
from numpy import linalg as LA
from numpy.linalg import matrix_power
import time 
import datetime 
import sys
import os
import random
import math
import scipy as sp
import scipy.linalg
from scipy.linalg import expm
from matplotlib.pyplot import *
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

startTime = time.time()
print ("STARTED:" , datetime.datetime.now().strftime("%d %B %Y %H:%M:%S"))

if len(sys.argv) < 6:
  print("Usage: python", str(sys.argv[0]), "READ-IN? " "SAVE-or-NOT? " "NCOL " "ITERS " "Symmetric case?")
  sys.exit(1)

READIN = int(sys.argv[1])
SAVE = int(sys.argv[2])
NCOL = int(sys.argv[3])
Niters_sim = int(sys.argv[4])
ZSYM = int(sys.argv[5])

if ZSYM == 1: 
    g, h = 1.0, 1.0
else: 
    g, h = 1/30., 1/15.
    
NMAT = 2
GAP = 1
dt = 1e-4
nsteps = int(1e-2/dt)
cut = int(0.25*Niters_sim) 
X = np.zeros((NMAT, NCOL, NCOL), dtype=complex)
mom_X = np.zeros((NMAT, NCOL, NCOL), dtype=complex)
f_X = np.zeros((NMAT, NCOL, NCOL), dtype=complex)
X_bak = np.zeros((NMAT, NCOL, NCOL), dtype=complex)
HAM, expDH, trX, trX2, trX4, evals, evals1 = [], [], [], [], [], [], []

if NMAT > 2:
    print ("Not supported yet. Edit the code or contact the author")
    sys.exit(1)
if READIN not in [0,1]:
    print ("Wrong input for READIN")
    sys.exit(1)
if SAVE not in [0,1]:
    print ("Wrong input for SAVE")
    sys.exit(1)
if ZSYM not in [0,1]:
    print ("Wrong input for symmetry choice")
    sys.exit(1)

print ("Hoppe type %2.0f-matrix model"%NMAT)
print ("NCOL = " "%3.0f " ","  " g = %2.4f , h = %2.4f " % (NCOL, g, h))
print ("--------------------------------------")

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

def hamil(X,mom_X):
    ham = potential(X)
    for j in range (NMAT):
        ham += 0.50 * np.trace(np.dot(mom_X[j],mom_X[j])).real
    return ham
         
def potential(X):
    s1, s2 = 0.0, 0.0
    for i in range (NMAT):

        if ZSYM == 1:
            s1 += 0.50 * np.trace(X[i] @ X[i]) 
        else:  
            s1 -= 0.50 * np.trace(X[i] @ X[i])

        s1 += (g/4.0)* np.trace(X[i] @ X[i] @ X[i] @ X[i])

    co = np.dot(X[0],X[1]) - np.dot(X[1],X[0])
    tr = np.trace(np.dot(co,co))
    s2 -= 0.50 * h *(tr.real)
    return ((s1+s2).real)*NCOL

def force(X):

    f_X = np.zeros((NMAT, NCOL, NCOL), dtype=complex)
    if ZSYM == 1:
        f_X[0] = X[0] + g*(X[0] @ X[0] @ X[0])
        f_X[1] = X[1] + g*(X[1] @ X[1] @ X[1])
    else: 
        f_X[0] = -X[0] + g*(X[0] @ X[0] @ X[0])
        f_X[1] = -X[1] + g*(X[1] @ X[1] @ X[1])

    f_X[0] -= h*dagger(comm(X[1], comm(X[0], X[1])))
    f_X[1] -= h*dagger(comm(X[0], comm(X[1], X[0])))

    for j in range(NMAT):
        if np.allclose(f_X[j], dagger(f_X[j])) == False:
            f_X[j] = makeH(f_X[j])

    return (f_X)*NCOL

def leapfrog(X,dt):

    mom_X = refresh_mom()
    ham_init = hamil(X,mom_X)

    for j in range(NMAT):
        X[j] += mom_X[j] * dt * 0.50

    for j in range(NMAT):

        for i in range(1, nsteps):
            f_X = force(X)
            mom_X[j] -= f_X[j]*dt
            X[j] += mom_X[j]*dt

    f_X = force(X)

    for j in range(NMAT):
        mom_X[j] -= f_X[j] * dt
        X[j] += mom_X[j] * dt * 0.50

    ham_final = hamil(X,mom_X)
    
    return X, ham_init, ham_final

def update(X,acc_count):

    X_bak = copy_fields(X)
    X, start, end = leapfrog(X, dt)
    change = end - start
    expDH.append(np.exp(-1.0*change))
    if np.exp(-change) < random.uniform(0,1):
        X = rejected_go_back_old_fields(X_bak)
        print(("REJECT: deltaH = " "%10.7f " " startH = " "%10.7f" " endH = " "%10.7f" % (change, start, end)))
    else:
        print(("ACCEPT: deltaH = " "%10.7f " "startH = " "%10.7f" " endH = " "%10.7f" % (change, start, end)))
        acc_count += 1 

    w, v = LA.eigh(X[0]) 
    w1, v1 = LA.eigh(X[1])
    evals.append(w/NCOL) # Normalize eigenvalues by N
    evals1.append(w1/NCOL)

    tmp0, tmp1 = np.trace(X[0]).real, np.trace(X[1]).real   
    trX.append((tmp0+tmp1)*0.5/NCOL)

    tmp2, tmp3 = np.trace(X[0] @ X[0]).real, np.trace(X[1] @ X[1]).real   
    trX2.append((tmp2+tmp3)*0.5/NCOL) 

    tmp4, tmp5 = np.trace(X[0] @ X[0] @ X[0] @ X[0]).real, np.trace(X[1] @ X[1] @ X[1] @ X[1]).real 
    trX4.append((tmp4+tmp5)*0.5/NCOL)

    if MDTU%GAP == 0:
        f3.write("%4.8f  \t %4.8f \n" %(tmp0/NCOL, tmp1/NCOL))
        f4.write("%4.8f  \t %4.8f \n" %(tmp2/NCOL, tmp3/NCOL))
        f5.write("%4.8f  \t %4.8f \n" %(tmp4/NCOL, tmp5/NCOL))
        
    return X,acc_count 


if __name__ == '__main__':
    
    if READIN == 0:
        print ("Starting from fresh")
        for i in range (NMAT):
            X[i] = 0.0  
                    
    if READIN == 1:
        name_f = "config_2MM_N{}.npy".format(NCOL)

        if os.path.isfile(name_f) == True: 
            print ("Reading old configuration file:", name_f)
            A = np.load(name_f)
            for i in range (NMAT):
                for j in range (NCOL):
                    for k in range (NCOL):
                        X[i][j][k] = A[i][j][k]

            for j in range(NMAT):
                if np.allclose(X[j], dagger(X[j])) == False:
                    print ("Input configuration 'X' not hermitian, ", LA.norm(X[j] - dagger(X[j])), "making it so")
                    X[j] = makeH(X[j])
        else:
            print ("Configuration not found, loading fresh")
            for i in range (NMAT):
                X[i] = 0.0

    f3 = open('t1_N%s_g%s_h%s.txt' %(NCOL,round(g,4),round(h,4)), 'w')
    f4 = open('t2_N%s_g%s_h%s.txt' %(NCOL,round(g,4),round(h,4)), 'w')
    f5 = open('t4_N%s_g%s_h%s.txt' %(NCOL,round(g,4),round(h,4)), 'w')

    acc_count = 0.0
    for MDTU in range (1, Niters_sim+1):

        X,acc_count = update(X,acc_count)

        if MDTU%10 == 0 and SAVE == 1:
            name_f = "config_2MM_N{}.npy".format(NCOL)
            print ("Saving configuration file: ", name_f)
            np.save(name_f, X)

    f3.close()
    f4.close()
    f5.close()
   
    if READIN == 0:
        expDH = expDH[cut:] 

    if acc_count/Niters_sim < 0.50:
        print("WARNING: Acceptance rate is below 50%")

    t2t4_plot = plt.figure(1) 
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.grid(which='major', axis='y', linestyle='--')
    MDTU = np.linspace(0, int(Niters_sim/GAP), int(Niters_sim/GAP), endpoint=True)
    plt.xlabel('Time units')

    if ZSYM == 1 and g == 1 and h == 1:
        plot(MDTU, trX2, 'blue', label=r'$t_{2}$')
        plot(MDTU, trX4, 'red', label=r'$t_{4}$')
        plt.axhline(y=0.421783612, color='blue', linestyle='--')
        plt.axhline(y=0.333341358, color='red', linestyle='--') 
        # Bootstrap results for reference

    else:
        plot(MDTU, trX2, 'blue', label=r'$t_{2}$')
        plot(MDTU, trX4, 'red', label=r'$t_{4}$')

    plt.legend(loc='best')
    evals = np.reshape(evals, (NCOL*Niters_sim)) 
    evals1 = np.reshape(evals1, (NCOL*Niters_sim))
    hist_plot = plt.figure(2)
    plt.hist(evals,color = "lightsalmon", density=False, bins=300) 
    plt.hist(evals1, color = '#6698EF', density=False, bins=300)  
    plt.ylabel(r'$\rho(x)$')
    plt.xlabel(r'$x$ (normalized by 1/N)')
    file="plots_2MM_N{}.pdf".format(NCOL)
    pp = PdfPages(file)
    pp.savefig(t2t4_plot, dpi = 300, transparent = True)
    pp.savefig(hist_plot, dpi = 300, transparent = True)
    pp.close()
    print("<Tr X / NCOL>", np.mean(trX), "+/-", (np.std(trX)/np.sqrt(np.size(trX) - 1.0)))
    print("<Tr X^2 / NCOL>", np.mean(trX2), "+/-", (np.std(trX2)/np.sqrt(np.size(trX2) - 1.0)))
    print("<Tr X^4 / NCOL>", np.mean(trX4), "+/-", (np.std(trX4)/np.sqrt(np.size(trX4) - 1.0)))
    print("<exp(-deltaH)>", np.mean(expDH), "+/-", np.std(expDH)/np.sqrt(np.size(expDH) - 1.0))
    print ("COMPLETED:" , datetime.datetime.now().strftime("%d %B %Y %H:%M:%S")) 
    endTime = time.time() 
    print ("Running time:", round(endTime - startTime, 2),  "seconds")