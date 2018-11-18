import random as r
import math as m
import numpy as np
x=[0.25,0.75,1.25,1.75,2.25]
y=[0.28,0.57,0.68,0.74,0.79]
a0=1.0
a1=1.0
while(True):
	b=[]
	c=[]
	f=[]
	d=[]
	for i in x:
		b.append(1-m.exp(-a1*i))
		c.append(a0*i*m.exp(-a1*i))
		f.append(a0*(1-m.exp(-a1*i)))
	z0=np.column_stack((b,c))	
	#print(z0)
	z0t=z0.transpose()
	e=np.dot(z0t,z0)
	#print(z0t)
	#print(e)
	I=np.linalg.inv(e)
	#print(i)
	for i in range(len(y)):
		d.append(y[i]-f[i])
	d=np.array(d)
	d=d.transpose()
	#print(z0t)
	#print(d)
	g=np.dot(z0t,d)
	dela=np.dot(I,g)
	sum=0
	for j in d:
		sum=sum+j*j
	#print(sum)
	if sum<0.0007:
		print("final fitted valuse for a0 and a1 are {:.6f} and {:.6f}".format(a0,a1))
		break
	else :
		a0=a0+dela[0]
		a1=a1+dela[1]			


#monte carlo using computer generated random numbers
#iterations=int(input("number of iterations: "))
iterations=1000
cir=0
sq=0
m=[]
n=[]
cirx=[]
ciry=[]
for i in range(iterations):
	x=r.random()
	y=r.random()
	m.append(x)
	n.append(y)
	a=(x*x)+(y*y)
	if a <= 1:
		cirx.append(x)
		ciry.append(y)
		cir+=1
	sq=i
pi=4*(cir/sq)
print("{:0.6f}".format(pi))


#monte carlo using manually generated random numbers (32 bit)
x0 = 12356789
y0 = 87238971
c = 16807
N = 2147483641 #($2^31 - 1)

for i in range(50000):
	x = ((c * x0)%N)
	y = ((c * y0)%N)
	x0 = x
	y0 = y
	x = x/N
	y = y/N
	#print("{:0.6f}".format(x)," {:0.6f}".format(y))
	a=(x*x)+(y*y)
	if a <= 1:
		#cirx.append(x)
		#ciry.append(y)
		cir+=1
	sq=i
pi=4*(cir/sq)
print("{:0.6f}".format(pi))