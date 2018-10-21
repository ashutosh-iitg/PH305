'''~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Lab No: 09
! Title : Difference & Shooting Method for Boundary Value Problems
! Date: 22/10/2018
! Name : Ashutosh Kumar Mandal
! Roll No: 160121010
! Webmail : ashutosh3004@iitg.ac.in
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''

import math
import copy
import random
import numpy as np
import matplotlib.pyplot as plt
from assignment4 import gaussSeidel

np.set_printoptions(formatter={'float': lambda x: "{0:0.6f}".format(x)})
degree_sign = u'\N{DEGREE SIGN}'

def act_f(x):
    return ((73.4523 * math.exp(0.1 * x)) - (53.4523 * math.exp(-0.1 * x)) + 20)

def grad_f(T,z,h_prime = 0.01,T_a = 20):
    grad_T = z
    grad_z = h_prime * (T - T_a)
    return grad_T,grad_z

def runge_kutta(T0,x0,z0,h = 0.5):
	lim = 10
	sol = []
	while(x0 <= lim):
		sol.append([x0,T0,z0])
		(fT,fz) = grad_f(T0,z0)
		kT1 = h * fT
		kz1 = h * fz
		(fT,fz) = grad_f(T0 + (0.5 * kT1),z0 + (0.5 * kz1))
		kT2 = h * fT
		kz2 = h * fz
		(fT,fz) = grad_f(T0 + (0.5 * kT2),z0 + (0.5 * kz2))
		kT3 = h * fT
		kz3 = h * fz
		(fT,fz) = grad_f(T0 + kT3,z0 + kz3)
		kT4 = h * fT
		kz4 = h * fz
		x1 = x0 + h
		T1 = T0 + (kT1 + (2 * kT2) + (2 * kT3) + kT4)/6
		z1 = z0 + (kz1 + (2 * kz2) + (2 * kz3) + kz4)/6
		x0 = x1
		T0 = T1
		z0 = z1
	return sol

def shooting(x1,x2,T1,T2):
	z0 = 10
	shot_1 = runge_kutta(T1,x1,z0)

	z0 = 20
	shot_2 = runge_kutta(T1,x1,z0)

	z0 = 10 + ((20 - 10)/(shot_2[len(shot_2) - 1][1] - shot_1[len(shot_1) - 1][1])) * (200 - shot_1[len(shot_1) - 1][1])
	hit = runge_kutta(T1,x1,z0)

	return shot_1,shot_2,hit

def difference(x1,x2,T1,T2, n = 10):
	h_prime = 0.01
	T_a = 20

	h = x2//n

	A = np.zeros((n-1,n-1))
	b = np.zeros(n-1)
	T = np.zeros(n+1)

	for i in range(n-1):
		A[i][i] = h**2 * h_prime + 2
		b[i] = h**2 * h_prime * T_a

		if (i != 0):
			A[i][i-1] = -1

		if (i != n-2):
			A[i][i+1] = -1

	b[0] += T1
	b[n-2] += T2
	#print(A,"\n\n",b)
	guess = random.sample(range(1,100),9)
	guess = np.array(guess)
	#print(guess)

	(sol,r) = gaussSeidel(A,b,guess)
	T[0] = T1
	T[n] = T2

	for i in range(n-1):
		T[i+1] = sol[i]

	#print(T)
	i = 0
	x = []
	while (i <= x2):
		x.append(i)
		i += h

	return T,x

if __name__ == "__main__" :

	(shot_1,shot_2,hit) = shooting(0,10,40,200)

	shot_1_x = []
	shot_1_T = []
	shot_2_x = []
	shot_2_T = []
	hit_x = []
	hit_T = []

	for j in shot_1:
		shot_1_x.append(j[0])
		shot_1_T.append(j[1])

	for j in shot_2:
		shot_2_x.append(j[0])
		shot_2_T.append(j[1])

	for j in hit:
		hit_x.append(j[0])
		hit_T.append(j[1])

	act_T = []
	t = 0
	h = 0.02
	act_x = []

	lim = 10
	while (t <= 10):
		act_T.append(act_f(t))
		act_x.append(t)
		t += h

	#Plot Solution of Shooting Method
	plt.grid(True)
	plt.plot(shot_1_x, shot_1_T, label = "Shot 1")
	plt.plot(shot_2_x, shot_2_T, label = "Shot 2")
	plt.plot(hit_x, hit_T, label = "Hit")
	plt.plot(act_x, act_T,'--', label = "Analytical Solution")
	plt.xlabel("Length (m)")
	plt.ylabel("Temperature ("+ degree_sign +"C)")
	plt.title("Shooting Method")
	plt.legend()
	plt.show()

	(diff_T, x) = difference(0,10,40,200)
	#print (diff_T)

	plt.grid(True)
	plt.plot(x,diff_T, label = "Diff. Method")
	plt.plot(act_x, act_T,'--', label = "Analytical Solution")
	plt.xlabel("Length (m)")
	plt.ylabel("Temperature ("+ degree_sign +"C)")
	plt.title("Finite Difference Method")
	plt.legend()
	plt.show()
	plt.show()