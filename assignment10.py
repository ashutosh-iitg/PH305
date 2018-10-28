import math
import matplotlib.pyplot as plt

def func(x,y):
	return eval("(2*x*y) + (2*x) - (x**2) - (2*y**2) + 72",{'x' : x, 'y' : y})

def trapezoid(a,b,c,d,N = 10):
	h = (b - a) / N
	k = (d - c) / N
	area = (b - a) * (d - c)
	T = 0
	for i in range (N + 1):
		yi = c + (i) * k
		for j in range (N + 1):
			xj = a + (j) * h
			#print (xj)
			if (j == 0 or j == N):
				if (i == 0 or i == N):
					T += func(xj,yi)
				else:
					T += 2 * func(xj,yi)
			elif (i == 0 or i == N):
				T += 2 * func(xj,yi)
			else:
				T += 4 * func(xj,yi)
			#print(T)
	T = (((h * k)/ 4) * T) / area
	return T


if __name__ == "__main__":
	print(" Average Temperature using Trapezoid method : {:0.6f} \n".format(trapezoid(0,8,0,6))) # a,b represents limits of x and c,d represents limits of y.
