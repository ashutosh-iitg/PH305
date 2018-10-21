'''=================================================
! Lab No: 04
! Title : Gauss Seidel and Jacobi Iteration method
! Date: 27/08/2018
! Name : Ashutosh Kumar Mandal
! Roll No: 160121010
! Webmail : ashutosh3004@iitg.ac.in
===================================================='''


import numpy as np
from matplotlib import pyplot as plt

def rowDominant(A,b):
	n = len(A)
	for i in range(n):
		maxindex = np.argmax(A[i])
		if ((np.sum(A[i]) - A[i,maxindex]) > A[i,maxindex]):
			raise valueError(" Matrix is not diagonal dominant. ")
		if (i != maxindex):
			A[[i,maxindex]] = A[[maxindex,i]]
			b[[i,maxindex]] = b[[maxindex,i]]
	return (A,b)


def gaussSeidel(A,b,guess):
	n = len(A)
	x = np.zeros(n)
	x_prev = np.zeros(n)
	iteration = 0
	error = 200
	x_gauss = []
	for i in range(n):
		x_prev[i] = guess[i]
	while( error > 0.0000000001 ):
		error = 0;
		iteration += 1
		i = 0
		for i in range(n):
			add = 0
			for j in range(n):
				if j == i:
					continue
				add = add + A[i][j] * x[j]
			x[i] = (b[i] - add)/A[i][i]
			error = abs(error + abs(x[i]) - abs(x_prev[i]))
			x_prev[i] = x[i]
		#print(x, error)
		x_gauss.append([iteration,x[0],x[1],x[2]])
		#print (iteration)

	return (x,x_gauss)


def jacobi(A,b,guess):
	n = len(A)
	x = np.zeros(n)
	iteration = 0
	error = 200
	x_jacobi = []
	while( error > 0.0000000001 ):
		iteration += 1
		i = 0
		for i in range(n):
			add = 0
			for j in range(n):
				if j == i:
					continue
				add = add + A[i][j] * guess[j]
			x[i] = (b[i] - add)/A[i][i]
		error = abs(np.sum(abs(x)) - np.sum(abs(guess)))
		for i in range(0,n):
			guess[i] = x[i]
		#print(x, error)
		x_jacobi.append([iteration,x[0],x[1],x[2]])
		#print (iteration)

	return (x,x_jacobi)



if __name__ == '__main__':
	#file = input(" Enter the file name: \n")
	data = np.genfromtxt("file04.csv",delimiter = ",")
	
	A = data[:,0:-1]
	b = data[:,-1]
	print(A,"\n\n")
	print(b,"\n\n")
	(A,b) = rowDominant(A,b)
	print(A,"\n\n")
	print(b,"\n\n")

	#guess = input(" Enter your guess values (in numpy array format) : \n")
	guess = np.array([1,1,1])

	(x1,x_gauss) = gaussSeidel(A, b,guess)
	(x2,x_jacobi) = jacobi(A, b, guess)

	print(" Solutions from Gauss-Seidel Iteration method : \n", x1);
	print(" Solutions from Jacobi Iteration method : \n", x2);
	#print(x_gauss[0])

	n_iters_g = []
	y1_g = []
	y2_g = []
	y3_g = []

	for i in x_gauss:
		n_iters_g.append(i[0])
		y1_g.append(i[1])
		y2_g.append(i[2])
		y3_g.append(i[3])

	n_iters_j = []
	y1_j = []
	y2_j = []
	y3_j = []

	for j in x_jacobi:
		n_iters_j.append(j[0])
		y1_j.append(j[1])
		y2_j.append(j[2])
		y3_j.append(j[3])

	plt.plot(n_iters_j, y1_j,"co-.")
	plt.plot(n_iters_j, y2_j,"yo-.")
	plt.plot(n_iters_j, y3_j,"ko-.")
	plt.plot(n_iters_g, y1_g,"bo-")
	plt.plot(n_iters_g, y2_g,"go-")
	plt.plot(n_iters_g, y3_g,"ro-")
	plt.show()