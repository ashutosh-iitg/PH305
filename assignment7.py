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


def act_f(x):
    return (4 * (math.exp(0.8 * x) - math.exp(-0.5 * x)) / 1.3) + 2 * math.exp(-0.5 * x)

def grad_f(x,y):
    return (4 * math.exp(0.8 * x) - (0.5 * y))

def euler(x0,y0,h):
    sol = []
    act_sol = []
    f =  open("euler.txt",'w')
    g =  open("actual.txt",'w')
    while(x0 <= 4):
        m = grad_f(x0,y0)
        y1 = y0 + m * h
        x1 = x0 + h
        sol.append([x1,y1])
        act_sol.append([x1,act_f(x1)])
        #print("{:.6f}".format(x1),"{:.6f}".format(y1),file = f)
        #print("{:.6f}".format(x1),"{:.6f}".format(act_f(x1)),file = g)
        x0 = x1
        y0 = y1
    f.close()
    g.close()
    return sol,act_sol

def midpoint(x0,y0,h):
    sol = []
    f =  open("midpoint.txt",'w')
    while (x0 <= 4):
        yh = y0 + grad_f(x0,y0) * h/2
        xh = x0 + h/2
        y1 = y0 + grad_f(xh,yh) * h
        x1 = x0 + h
        sol.append([x1,y1])
        #print("{:.6f}".format(x1),"{:.6f}".format(y1),file = f)
        x0 = x1
        y0 = y1
    f.close()
    return sol

def heun(x0,y0,h):
    sol = []
    f = open("heun.txt",'w')
    while (x0 <= 4):
        yh = y0 + grad_f(x0,y0) * h
        xh = x0 + h
        y1 = y0 + (grad_f(x0,y0) + grad_f(xh,yh)) * h / 2
        x1 = x0 + h
        sol.append([x1,y1])
        #print("{:.6f}".format(x1),"{:.6f}".format(y1),file = f)
        x0 = x1
        y0 = y1
    f.close()
    return sol



if __name__ == "__main__" :
    
    (sol_euler,act_sol) = euler(0,2,0.1)
    sol_midpoint = midpoint(0,2,0.1)
    sol_heun = heun(0,2,0.1)
    
    x_act = []
    y_act = []
    x_euler = []
    y_euler = []
    x_midpoint = []
    y_midpoint = []
    x_heun = []
    y_heun = []
    
    for j in act_sol:
        x_act.append(j[0])
        y_act.append(j[1])
    
    for j in sol_euler:
        x_euler.append(j[0])
        y_euler.append(j[1])
    
    for j in sol_midpoint:
        x_midpoint.append(j[0])
        y_midpoint.append(j[1])
        
    for j in sol_heun:
        x_heun.append(j[0])
        y_heun.append(j[1])
    
    plt.plot(x_euler,y_euler,"m-",label = "Euler Approximation")
    plt.plot(x_midpoint,y_midpoint,"r-.",label = "Mid-point Approximation")
    plt.plot(x_heun,y_heun,"g--",label = "Heun Approximation")
    plt.plot(x_act,y_act,"b--",label = "Actual Solution")
    plt.xlabel("x -->")
    plt.ylabel("y -->")
    plt.legend()
    plt.show()
