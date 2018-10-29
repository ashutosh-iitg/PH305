'''~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Lab No: 10
! Title : Numerical methods for definite inegrals
! Date: 29/10/2018
! Name : Ashutosh Kumar Mandal
! Roll No: 160121010
! Webmail : ashutosh3004@iitg.ac.in
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''

import math
import numpy as np

np.set_printoptions(formatter={'float': lambda x: "{0:0.2f}".format(x)})

def func(x,y):
	return eval("(2*x*y) + (2*x) - (x**2) - (2*y**2) + 72",{'x' : x, 'y' : y})

def trapezoid(a,b,c,d,N = 10):
	h = (b - a) / N
	k = (d - c) / N
	area = (b - a) * (d - c)
	T = 0
	for j in range (N + 1):
		yj = c + (j * k)
		for i in range (N + 1):
			xi = a + (i * h)
			#print (xi)
			if (j == 0 or j == N):
				if (i == 0 or i == N):
					T += func(xi,yj)
				else:
					T += 2 * func(xi,yj)
			elif (i == 0 or i == N):
				T += 2 * func(xi,yj)
			else:
				T += 4 * func(xi,yj)
			#print(T)
	T = (((h * k)/ 4) * T) / area
	return T

def simpson_onethird(a,b,c,d,N = 2):
	area = (b - a) * (d - c)
	h = (b - a) / N
	k = (d - c) / N

	S = np.zeros((N + 1,N + 1))
	F = np.zeros((N + 1,N + 1))
	T = 0

	S[0][0] = 1
	S[0][N] = 1
	S[N][0] = 1
	S[N][N] = 1
	for i in range (1,N) : 
		if (i % 2) == 0:
			S[0][i] = 2
			S[i][0] = 2
			S[N][i] = 2
			S[i][N] = 2
		else:
			S[0][i] = 4
			S[i][0] = 4
			S[N][i] = 4
			S[i][N] = 4
	for i in range(1,N):
		for j in range (1,N):
			S[i][j] = S[i][0] * S[0][j]

	for i in range(N + 1):
		y = c + i * k
		for j in range (N + 1):
			x = a + j * h
			F[i][j] = func(x,y)
	
	for i in range(N + 1):
		for j in range (N + 1):
			T += S[i][j] * F[i][j]

	T = (T * h * k) / (9 * area)
	return T

def simpson_three_eighth(a,b,c,d,N = 3):
	area = (b - a) * (d - c)
	h = (b - a) / N
	k = (d - c) / N

	S = np.zeros((N + 1,N + 1))
	F = np.zeros((N + 1,N + 1))
	T = 0

	S[0][0] = 1
	S[0][N] = 1
	S[N][0] = 1
	S[N][N] = 1
	for i in range (1,N) : 
		if (i % 3) == 0:
			S[0][i] = 2
			S[i][0] = 2
			S[N][i] = 2
			S[i][N] = 2
		else:
			S[0][i] = 3
			S[i][0] = 3
			S[N][i] = 3
			S[i][N] = 3

	for i in range(1,N):
		for j in range (1,N):
			S[i][j] = S[i][0] * S[0][j]

	for i in range(N + 1):
		y = c + i * k
		for j in range (N + 1):
			x = a + j * h
			F[i][j] = func(x,y)
	
	for i in range(N + 1):
		for j in range (N + 1):
			T += S[i][j] * F[i][j]

	T = (h * k * T * 9) / (64 * area)
	return T

if __name__ == "__main__":
	N = int(input("Enter no. of subdivisions of the domain : "))
	print(" Average Temperature using Trapezoid method : {:0.6f} \n".format(trapezoid(0,8,0,6,N))) # a,b represents limits of x and c,d represents limits of y.
	
	while (N%2 != 0):
		N = int(input("Enter no. of subdivisions of the domain : "))
	print(" Average Temperature using Simpson 1/3 method : {:0.6f} \n".format(simpson_onethird(0,8,0,6,N)))

	while (N%3 != 0):
		N = int(input("Enter no. of subdivisions of the domain : "))
	print(" Average Temperature using Simpson 3/8 method : {:0.6f} \n".format(simpson_three_eighth(0,8,0,6)))