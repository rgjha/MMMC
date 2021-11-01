#!/usr/bin/python3
# -*- coding: utf-8 -*-
import time 
import datetime 
import sys
import numpy as np
import random
import math
import os 
from numpy import linalg as LA
from matplotlib.pyplot import *
from matplotlib import pyplot as plt

startTime = time.time()
print ("STARTED:" , datetime.datetime.now().strftime("%d %B %Y, %H:%M:%S"))
if len(sys.argv) < 5:
  print("Usage:python",str(sys.argv[0]),"READ-IN? " "SAVE-or-NOT? " "NCOL " "NITERS")
  sys.exit(1)

READIN = int(sys.argv[1])
SAVE = int(sys.argv[2])
NCOL = int(sys.argv[3]) 
Niters_sim = int(sys.argv[4])

NMAT = 3
g = 1.0
dt = 1e-3
nsteps = int(0.5/dt) 
GAP = 2.
cut = int(0.25*Niters_sim) 

if Niters_sim%GAP != 0:
  print("'Niters_sim' mod 'GAP' is not zero ")
  sys.exit(1) 
if READIN not in [0,1]:
    print ("Wrong input for READIN")
    sys.exit(1)
if SAVE not in [0,1]:
    print ("Wrong input for SAVE")
    sys.exit(1)
  
X = np.zeros((NMAT, NCOL, NCOL), dtype=complex)
mom_X = np.zeros((NMAT, NCOL, NCOL), dtype=complex)
f_X = np.zeros((NMAT, NCOL, NCOL), dtype=complex)
X_bak = np.zeros((NMAT, NCOL, NCOL), dtype=complex)
HAM, expDH, trX2, trX4, MOM = [], [], [], [], []
t2_ex = [None] * NMAT
t4_ex = [None] * NMAT

print ("Matrix integral simulation of%2.0f MM"%(NMAT)) 
print ("NCOL =" " %3.0f " ","  " and g =" " %4.2f" % (NCOL, g)) 
print ("------------------------------------------------------")

def dagger(a):
	return np.transpose(a).conj()

def box_muller():  
	PI = 2.0*math.asin(1.0);    
	r = random.uniform(0,1)
	s = random.uniform(0,1)
	p = np.sqrt(-2.0*np.log(r)) * math.sin(2.0*PI*s)
	q = np.sqrt(-2.0*np.log(r)) * math.cos(2.0*PI*s)
	return p,q

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
		print ("WARNING: Couldn't make hermitian")
	return tmp2

def hamil(X,mom_X):
	ham = potential(X) 
	for j in range (NMAT):
		ham += 0.50 * np.trace(np.dot(mom_X[j],mom_X[j])).real 
	return ham  

def potential(X):
	pot = 0.0 
	for i in range (NMAT):
		pot += 0.50 * np.trace(np.dot(X[i],X[i])).real   
		pot += (g/4.0)* np.trace(X[i] @ X[i] @ X[i] @ X[i]).real
	return pot*NCOL

def force(X): 
	for i in range (NMAT): 
		f_X[i] = (X[i] + (g*(X[i] @ X[i] @ X[i])))*NCOL
	for j in range(NMAT):
		if np.allclose(f_X[j], dagger(f_X[j])) == False:
			f_X[j] = makeH(f_X[j])
	return f_X

def leapfrog(X,dt):

	mom_X = refresh_mom()
	ham_init = hamil(X,mom_X)

	for j in range(NMAT):
		X[j] += mom_X[j] * dt * 0.5 

	for i in range(1, nsteps+1):
		f_X = force(X)
		for j in range(NMAT):
			mom_X[j] -= f_X[j] * dt
			X[j] += mom_X[j] * dt

	f_X = force(X)
	for j in range(NMAT):
		mom_X[j] -= f_X[j] * dt
		X[j] += mom_X[j] * dt * 0.5 

	ham_final = hamil(X,mom_X)
	return X, ham_init, ham_final

def update(X, acc_count):
	
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

	if MDTU%GAP == 0:
		t2_ex[0] = np.trace(np.dot(X[0],X[0])).real
		trX2.append(t2_ex[0]/NCOL)
		t4_ex[0] = np.trace((X[0] @ X[0] @ X[0] @ X[0])).real
		trX4.append(t4_ex[0]/NCOL)
		
		if NMAT > 1:
			for i in range (1, NMAT):
				t2_ex[i] = np.trace(np.dot(X[i],X[i])).real
				t4_ex[i] = np.trace((X[i] @ X[i] @ X[i] @ X[i])).real
		
		for item in t2_ex:
			f3.write("%4.8f " % (item/NCOL))
		for item in t4_ex:
			f4.write("%4.8f " % (item/NCOL))
		f3.write("\n")
		f4.write("\n")

	return X, acc_count

if __name__ == '__main__':

	if READIN == 0:
		print ("Loading fresh configuration")
		for i in range (NMAT): 
			X[i] = 0.0  

	if READIN == 1:

		name_f = "config_1MM_N{}.npy".format(NCOL)

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
			for i in range (NMAT): 
				print ("Configuration not found, loaded fresh")
				X[i] = 0.0


	f3 = open('t2_1MM_N%s_g%s.txt' %(NCOL,round(g,4)), 'w')
	f4 = open('t4_1MM_N%s_g%s.txt' %(NCOL,round(g,4)), 'w')

	acc_count = 0.
	for MDTU in range (1, Niters_sim+1):

		X, acc_count = update(X, acc_count)
		if MDTU%10 == 0 and SAVE == 1:
			name_f = "config_1MM_N{}.npy".format(NCOL)
			print ("Saving configuration file: ", name_f)
			np.save(name_f, X)

	f3.close()
	f4.close()

	if NMAT == 1:  
		t2_exact = (((12*g)+1)**(3/2.) - 18*g - 1)/(54*g*g)
		# Exact result for 1MM quartic potential with g = 1
		plt.rc('text', usetex=True)
		plt.rc('font', family='serif')
		MDTU = np.linspace(0, int(Niters_sim/GAP), int(Niters_sim/GAP), endpoint=True)
		plt.ylabel(r'Tr(X$^2)/N$',fontsize=12)
		plt.xlabel('Time units', fontsize=12)
		plt.grid(which='major', axis='y', linestyle='--')
		plt.axhline(y=t2_exact, color='teal', linestyle='--')
		plt.figure(1)
		plot (MDTU, trX2, 'teal') 
		outname = "1MM_N%s_g%s" %(NCOL, g)
		plt.savefig(outname+'.pdf')
	print ("------------------------------------------------------")
	print ("Acceptance rate: ", (acc_count/Niters_sim)*100,"%") 
	if acc_count/Niters_sim < 0.50:
		print("WARNING: Acceptance rate is below 50%")

	if READIN == 0:
		trX2 = trX2[cut:]
		trX4 = trX4[cut:]
		expDH = expDH[cut:] 

	print ("COMPLETED:" , datetime.datetime.now().strftime("%d %B %Y, %H:%M:%S"))
	endTime = time.time() 
	print ("Running time:", round(endTime - startTime, 2),  "seconds")