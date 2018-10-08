'''~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Lab No: 07
! Title : Euler and RK-4 Method for Harmonic Oscillator
! Date: 08/10/2018
! Name : Ashutosh Kumar Mandal
! Roll No: 160121010
! Webmail : ashutosh3004@iitg.ac.in
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''

import math
import copy
import numpy as np
import matplotlib.pyplot as plt

def act_f(t):
    return 1 * math.cos(t)

def grad_f(x,u):
    grad_x = u
    grad_u = -x
    return grad_x,grad_u

def oscillator(x,u,m = 1,k = 1):
	p = m * u
	E = (p**2 / (2 * m)) + (k * x**2) / 2
	return p,E

def euler(t0,x0,u0,h = 0.02 * 2 * math.pi):
    lim = 6 * math.pi
    sol = []
    act_sol = []
    f =  open("euler.txt",'w')
    g =  open("actual.txt",'w')
    while(t0 <= lim):
        (mx,mu) = grad_f(x0,u0)
        u1 = u0 + mu * h
        x1 = x0 + mx * h
        t1 = t0 + h
        p,E = oscillator(x1,u1)
        sol.append([t1,x1,u1,p,E])
        act_sol.append([t1,act_f(t1)])
        #print("{:.6f}".format(p),"{:.6f}".format(x1),file = f)
        #print("{:.6f}".format(t1),"{:.6f}".format(act_f(t1)),file = g)
        t0 = t1
        x0 = x1
        u0 = u1
    f.close()
    g.close()
    return sol,act_sol

def runge_kutta(t0,x0,u0,h = 0.02 * 2 * math.pi):
	lim = 6 * math.pi
	sol = []
	f = open("runge_kutta.txt",'w')
	while(t0 <= lim):
		(fx,fu) = grad_f(x0,u0)
		kx1 = h * fx
		ku1 = h * fu
		(fx,fu) = grad_f(x0 + (0.5 * kx1),u0 + (0.5 * ku1))
		kx2 = h * fx
		ku2 = h * fu
		(fx,fu) = grad_f(x0 + (0.5 * kx2),u0 + (0.5 * ku2))
		kx3 = h * fx
		ku3 = h * fu
		(fx,fu) = grad_f(x0 + kx3,u0 + ku3)
		kx4 = h * fx
		ku4 = h * fu
		t1 = t0 + h
		x1 = x0 + (kx1 + (2 * kx2) + (2 * kx3) + kx4) / 6
		u1 = u0 + (ku1 + (2 * ku2) + (2 * ku3) + ku4) / 6
		p,E = oscillator(x1,u1)
		sol.append([t1,x1,u1,p,E])
		#print("{:.6f}".format(p),"{:.6f}".format(x1),file = f)
		t0 = t1
		x0 = x1
		u0 = u1
	f.close()
	return sol


if __name__ == "__main__" :
    
    (sol_euler,act_sol) = euler(0,1,0)
    sol_runge = runge_kutta(0,1,0)
    
    t_act = []
    x_act = []
    t_euler = []
    x_euler = []
    u_euler = []
    p_euler = []
    E_euler = []
    t_runge = []
    x_runge = []
    u_runge = []
    p_runge = []
    E_runge = []
    
    for j in act_sol:
        t_act.append(j[0])
        x_act.append(j[1])
    
    for j in sol_euler:
        t_euler.append(j[0])
        x_euler.append(j[1])
        u_euler.append(j[2])
        p_euler.append(j[3])
        E_euler.append(j[4])

    for j in sol_runge:
        t_runge.append(j[0])
        x_runge.append(j[1])
        u_runge.append(j[2])
        p_runge.append(j[3])
        E_runge.append(j[4])
    
    plt.rcParams["figure.figsize"] = [6,5.5]

    plt.plot(t_euler,x_euler,"m-",label = "Euler Approximation")
    plt.plot(t_runge,x_runge,"r.",label = "RK-4 Approximation")
    plt.plot(t_act,x_act,"g-",label = "Actual Solution")
    plt.xlabel("Time")
    plt.ylabel("Position")
    plt.legend()
    plt.grid(True)
    plt.title("Position vs Time for a Harmonic Oscillator")
    plt.show()

    plt.plot(x_euler,p_euler,"m-",label = "Euler Approximation")
    plt.plot(x_runge,p_runge,"b-.", label = "RK-4 Approximation")
    plt.xlabel("Position")
    plt.ylabel("Momentum")
    plt.legend()
    plt.grid(True)
    plt.title("Phase-Space Diagram for Harmonic Oscillator")
    plt.ylim(-3,3)
    plt.show()

    t_euler = np.array(t_euler) / ( 2 * math.pi)
    t_runge = np.array(t_runge) / ( 2 * math.pi)
    plt.plot(t_euler,E_euler,"m-",label = "Euler Approximation")
    plt.plot(t_runge,E_runge,"b-.",label = "RK-4 Approximation")
    plt.xlabel("t/T (T = 2$\pi$)")
    plt.ylabel("Energy")
    plt.legend()
    plt.grid(True)
    plt.title("Energy vs Time(t/T) for Harmonic Oscillator")
    plt.xlim(0,3.1)
    plt.show()