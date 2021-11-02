#!/usr/bin/python3
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
from matplotlib import pyplot

startTime = time.time()
print ("STARTED:" , datetime.datetime.now().strftime("%d %B %Y %H:%M:%S"))

if len(sys.argv) < 8:
  print("Usage: python", str(sys.argv[0]), "READ-IN? " "SAVE-or-NOT? " "NCOL " "NITER " "g "  "c " "kappa ")
  sys.exit(1)

READIN = int(sys.argv[1])
SAVE = int(sys.argv[2])
NCOL = int(sys.argv[3])
Niters_sim = int(sys.argv[4])
g = float(sys.argv[5])
c = float(sys.argv[6])
kappa = float(sys.argv[7])

if READIN not in [0,1]:
    print ("Wrong input for READIN")
    sys.exit(1)
if SAVE not in [0,1]:
    print ("Wrong input for SAVE")
    sys.exit(1)

NMAT = 4
dt = 1e-4   
nsteps = int(1e-2/dt)
GAP = 1

X = np.zeros((NMAT, NCOL, NCOL), dtype=complex)
mom_X = np.zeros((NMAT, NCOL, NCOL), dtype=complex)
f_X = np.zeros((NMAT, NCOL, NCOL), dtype=complex)
X_bak = np.zeros((NMAT, NCOL, NCOL), dtype=complex)
expDH = [] 

if NMAT == 3:
    trXsq, trYsq, trZsq = [], [], [] 
if NMAT == 4:
    trXsq, trYsq, trZsq, trWsq = [], [], [], [] 
if NMAT > 4:
    print ("This is not supported & checked yet! Lower NMAT or modify the code (not difficult!) or contact the author")
    sys.exit(1)


print ("Matrix chain with %2.0f matrices" %(NMAT))
print ("NCOL = " "%3.0f "","" and g = %4.2f,c = %4.2f,kappa = %4.2f"%(NCOL, g, c, kappa))
print ("----------------------------------------------------")

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
    s1 = 0.0
    for i in range (NMAT):
        s1 += 0.50 * np.trace(np.dot(X[i],X[i]))
        s1 += (g/4.0)* np.trace((matrix_power(X[i], 4)))
        if i == NMAT-1:
            pre = kappa
        else:
            pre = c
        s1 -= pre*np.trace(np.dot(X[i],X[(i+1)%NMAT])) 
    return (s1*NCOL).real 

def force(X):
    f_X = np.zeros((NMAT, NCOL, NCOL), dtype=complex)

    for i in range (NMAT):
        if i == NMAT-1:
            pre1 = kappa
            pre2 = c
        elif i == 0:
            pre1 = c
            pre2 = kappa
        else:
            pre1 = pre2 = c
        f_X[i] = (X[i]+(g*(matrix_power(X[i], 3))) - (pre1*X[(i+1)%NMAT]) - (pre2*X[(i-1+NMAT)%NMAT]))*NCOL

    for j in range(NMAT):
        if np.allclose(f_X[j], dagger(f_X[j])) == False:
            f_X[j] = makeH(f_X[j])

    return f_X 

def leapfrog(X,dt):
    mom_X = refresh_mom()
    ham_init = hamil(X,mom_X)

    for j in range(NMAT):
        X[j] = X[j] + (mom_X[j] * dt * 0.50)

    for j in range(NMAT):

        for i in range(1, nsteps):
            f_X = force(X)
            mom_X[j] = mom_X[j] - (f_X[j]*dt)
            X[j] = X[j] + (mom_X[j]*dt)

    f_X = force(X)

    for j in range(NMAT):

        mom_X[j] = mom_X[j] - (f_X[j] * dt)
        X[j] = X[j] + (mom_X[j] * dt * 0.50)

    ham_final = hamil(X,mom_X)
    
    return X, ham_init, ham_final

def update(X, acc_count):

    X_bak = copy_fields(X)
    X, start, end = leapfrog(X, dt)
    change = end - start
    if np.exp(-change) < random.uniform(0,1):
        X = rejected_go_back_old_fields(X_bak)
        print(("REJECT: deltaH = " "%10.7f " " startH = " "%10.7f" " endH = " "%10.7f" % (change, start, end)))
    else:
        print(("ACCEPT: deltaH = " "%10.7f " "startH = " "%10.7f" " endH = " "%10.7f" % (change, start, end)))
        acc_count += 1 

    if NMAT >= 3:
        tmp0 = np.trace(np.dot(X[0],X[0])).real
        tmp00 = np.trace(np.dot(np.dot(X[0],X[0]),np.dot(X[0],X[0]))).real
        trXsq.append(tmp0/NCOL)    
        tmp1 = np.trace(np.dot(X[1],X[1])).real
        tmp11 = np.trace(np.dot(np.dot(X[1],X[1]),np.dot(X[1],X[1]))).real
        trYsq.append(tmp1/NCOL)
        tmp2 = np.trace(np.dot(X[2],X[2])).real
        tmp22 = np.trace(np.dot(np.dot(X[2],X[2]),np.dot(X[2],X[2]))).real
        trZsq.append(tmp2/NCOL)

        if NMAT == 4:
            tmp3 = np.trace(np.dot(X[3],X[3])).real
            tmp33 = np.trace(np.dot(np.dot(X[3],X[3]),np.dot(X[3],X[3]))).real
            trWsq.append(tmp3/NCOL)

    if MDTU%GAP == 0:

        if NMAT == 3: 
            f3.write("%4.8f \t %4.8f \t %4.8f \n" %(tmp0/NCOL, tmp1/NCOL, tmp2/NCOL))
            f4.write("%4.8f \t %4.8f \t %4.8f \n" %(tmp00/NCOL, tmp11/NCOL, tmp22/NCOL))
        if NMAT == 4:
            f3.write("%4.8f \t %4.8f \t %4.8f \t %4.8f \n" %(tmp0/NCOL, tmp1/NCOL, tmp2/NCOL, tmp3/NCOL))
            f4.write("%4.8f \t %4.8f \t %4.8f \t %4.8f \n" %(tmp00/NCOL, tmp11/NCOL, tmp22/NCOL, tmp33/NCOL))

    return X, acc_count 


if __name__ == '__main__':
    
    
    if READIN == 0:
        print ("Loading fresh configuration")
        for i in range (NMAT):
            X[i] = 0.0

  
    if READIN == 1:
        name_f = "config_{}MM_N{}_g{}.npy".format(NMAT, NCOL, g)
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
            print ("Configuration not found, loaded fresh")
            for i in range (NMAT): 
                X[i] = 0.0


    f3 = open('t2_N%s_D%s_g%s.txt' %(NCOL,kappa,g), 'w')
    f4 = open('t4_N%s_D%s_g%s.txt' %(NCOL,kappa,g), 'w')

    acc_count = 0
    for MDTU in range (1, Niters_sim+1):
        X,acc_count = update(X,acc_count)

        if MDTU%10 == 0 and SAVE == 1:
            name_f = "config_{}MM_N{}_g{}.npy".format(NMAT, NCOL, g)
            print ("Saving configuration file: ", name_f)
            np.save(name_f, X)

    f3.close()
    f4.close()

    if acc_count/Niters_sim < 0.50:
        print("WARNING: Acceptance rate is below 50%")

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.grid(which='major', axis='y', linestyle='--')
    MDTU = np.linspace(0, int(Niters_sim/GAP), int(Niters_sim/GAP), endpoint=True)
    plt.ylabel(r'$\langle \rm{Tr}(X_{1,2,3}^{2})/N \rangle$')
    plt.xlabel('Time units')

    if NMAT == 3: 
        plot(MDTU, trXsq, 'red')
        plot(MDTU, trYsq, 'teal')
        plot(MDTU, trZsq, 'blue')

    if NMAT == 4:
        plot(MDTU, trXsq, 'red')
        plot(MDTU, trYsq, 'teal')
        plot(MDTU, trZsq, 'blue')
        plot(MDTU, trWsq, 'green')

    outname = "%sMM_N%s_k%s" %(NMAT, NCOL, kappa)
    plt.savefig(outname+'.pdf')
    print ("Fraction % of MDTU accepted = ", (acc_count/Niters_sim)*100) 
    print ("COMPLETED:" , datetime.datetime.now().strftime("%d %B %Y %H:%M:%S")) 
    endTime = time.time() 
    print ("Running time:", round(endTime - startTime, 2),  "seconds")