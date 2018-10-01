'''~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Lab No: 07
! Title : Euler, Mid-point and Heun's method
! Date: 01/10/2018
! Name : Ashutosh Kumar Mandal
! Roll No: 160121010
! Webmail : ashutosh3004@iitg.ac.in
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''

import math
import copy
from matplotlib import pyplot as plt

def act_f(t):
    return 1 * math.cos(t)

def grad_f(x,u):
    grad_x = u
    grad_u = -x
    return grad_x,grad_u

def euler(t0,x0,u0,h = 0.01):
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
        sol.append([t1,x1,u1])
        act_sol.append([t1,act_f(t1)])
        #print("{:.6f}".format(t1),"{:.6f}".format(x1),file = f)
        #print("{:.6f}".format(t1),"{:.6f}".format(act_f(t1)),file = g)
        t0 = t1
        x0 = x1
        u0 = u1
    f.close()
    g.close()
    return sol,act_sol

def midpoint(t0,x0,u0,h = 0.01):
    lim = 6 * math.pi
    sol = []
    f = open("midpoint.txt",'w')
    while (t0 <= lim):
        (mx,mu) = grad_f(x0,u0)
        uh = u0 + mu * h/2
        xh = x0 + mx * h/2
        (mxh,muh) = grad_f(xh,uh)
        u1 = u0 + muh * h
        x1 = x0 + mxh * h
        t1 = t0 + h
        sol.append([t1,x1,u1])
        #print("{:.6f}".format(t1),"{:.6f}".format(x1),file = f)
        t0 = t1
        x0 = x1
        u0 = u1
    f.close()
    return sol

def heun(t0,x0,u0,h = 0.01):
    lim = 6 * math.pi
    sol = []
    f = open("heun.txt",'w')
    while (t0 <= lim):
        (mx,mu) = grad_f(x0,u0)
        uh = u0 + mu * h
        xh = x0 + mx * h
        (mxh,muh) = grad_f(xh,uh)
        u1 = u0 + (muh + mu) * h / 2
        x1 = x0 + (mxh + mx) * h / 2
        t1 = t0 + h
        sol.append([t1,x1,u1])
        #print("{:.6f}".format(t1),"{:.6f}".format(x1),file = f)
        t0 = t1
        x0 = x1
        u0 = u1
    f.close()
    return sol


if __name__ == "__main__" :
    
    (sol_euler,act_sol) = euler(0,1,0)
    sol_midpoint = midpoint(0,1,0)
    sol_heun = heun(0,1,0)
    
    t_act = []
    x_act = []
    t_euler = []
    x_euler = []
    u_euler = []
    t_midpoint = []
    x_midpoint = []
    u_midpoint = []
    t_heun = []
    x_heun = []
    y_heun = []
    
    for j in act_sol:
        t_act.append(j[0])
        x_act.append(j[1])
    
    for j in sol_euler:
        t_euler.append(j[0])
        x_euler.append(j[1])
        u_euler.append(j[2])
    
    for j in sol_midpoint:
        t_midpoint.append(j[0])
        x_midpoint.append(j[1])
        u_midpoint.append(j[2])
        
    for j in sol_heun:
        t_heun.append(j[0])
        x_heun.append(j[1])
        y_heun.append(j[2])
    
    plt.plot(t_euler,x_euler,"m-",label = "Euler Approximation")
    #plt.plot(t_euler,u_euler)
    plt.plot(t_midpoint,x_midpoint,"b-.",label = "Mid-point Approximation")
    plt.plot(t_heun,x_heun,"g-.",label = "Heun Approximation")
    plt.plot(t_act,x_act,"r--",label = "Actual Solution")
    plt.xlabel("Time")
    plt.ylabel("Position")
    plt.legend()
    plt.show()
