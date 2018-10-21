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
import numpy as np
import matplotlib.pyplot as plt

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

def difference():
	return 0

def shooting(x1,x2,T1,T2):
	z0 = 10
	shot_1 = runge_kutta(T1,x1,z0)

	z0 = 20
	shot_2 = runge_kutta(T1,x1,z0)

	z0 = 10 + ((20 - 10)/(shot_2[len(shot_2) - 1][1] - shot_1[len(shot_1) - 1][1])) * (200 - shot_1[len(shot_1) - 1][1])
	hit = runge_kutta(T1,x1,z0)

	return shot_1,shot_2,hit

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
	x = []

	lim = 10
	while (t <= 10):
		act_T.append(act_f(t))
		x.append(t)
		t += h

	#Plot Solution of Shooting Method
	plt.grid(True)
	plt.plot(shot_1_x, shot_1_T, label = "Shot 1")
	plt.plot(shot_2_x, shot_2_T, label = "Shot 2")
	plt.plot(hit_x, hit_T, label = "Hit")
	plt.plot(x, act_T,'--', label = "Analytical Solution")
	plt.xlabel("Length (m)")
	plt.ylabel("Temperature ("+ degree_sign +"C)")
	plt.title("Shooting Method")
	plt.legend()
	plt.show()