import math
import copy
from matplotlib import pyplot as plt

def grad_f(y):
    return (-2.2067 * 1e-12 * (y ** 4 - (81 * 10 ** 8)))

def euler(x0,y0,h = 10):
    sol = []
    sol.append([x0,y0])
    while(x0 < 480):
        m = grad_f(y0)
        y1 = y0 + m * h
        x1 = x0 + h
        sol.append([x1,y1])
        x0 = x1
        y0 = y1
    print("Temperature at t = {:.1f} (Euler's method) : {:.6f}".format(x1,y1))
    return sol

def midpoint(x0,y0,h = 10):
    sol = []
    sol.append([x0,y0])
    while (x0 < 480):
        yh = y0 + grad_f(y0) * h/2
        xh = x0 + h/2
        y1 = y0 + grad_f(yh) * h
        x1 = x0 + h
        sol.append([x1,y1])
        x0 = x1
        y0 = y1
    print("Temperature at t = {:.1f} (Mid-Point method) : {:.6f}".format(x1,y1))
    return sol

def heun(x0,y0,h = 10):
    sol = []
    sol.append([x0,y0])
    while (x0 < 480):
        yh = y0 + grad_f(y0) * h
        xh = x0 + h
        y1 = y0 + (grad_f(y0) + grad_f(yh)) * h / 2
        x1 = x0 + h
        sol.append([x1,y1])
        x0 = x1
        y0 = y1
    print("Temperature at t = {:.1f} (Heun's method) : {:.6f}".format(x1,y1))
    return sol

def ralston(x0,y0,h = 10):
    sol = []
    sol.append([x0,y0])
    a1 = 1/3
    a2 = 2/3
    p1 = 3/4
    q1 = 3/4
    while (x0 < 480):
        k1 = grad_f(y0)
        k2 = grad_f(y0 + q1 * k1 * h)
        y1 = y0 + (a1 * k1 + a2 * k2) * h
        x1 = x0 + h
        sol.append([x1,y1])
        x0 = x1
        y0 = y1
    print("Temperature at t = {:.1f} (Ralston's method) : {:.6f}".format(x1,y1))
    return sol


if __name__ == "__main__" :

    h = float(input("Enter step size : \n "))
    
    sol_euler = euler(0,1200,h)
    sol_midpoint = midpoint(0,1200,h)
    sol_heun = heun(0,1200,h)
    sol_ralston = ralston(0,1200,h)
    
    x_euler = []
    y_euler = []
    x_midpoint = []
    y_midpoint = []
    x_heun = []
    y_heun = []
    x_ralston = []
    y_ralston = []
    
   
    for j in sol_euler:
        x_euler.append(j[0])
        y_euler.append(j[1])
    
    for j in sol_midpoint:
        x_midpoint.append(j[0])
        y_midpoint.append(j[1])
        
    for j in sol_heun:
        x_heun.append(j[0])
        y_heun.append(j[1])

    for j in sol_ralston:
        x_ralston.append(j[0])
        y_ralston.append(j[1])
    
    plt.plot(x_euler,y_euler,"m-",label = "Euler Approximation")
    plt.plot(x_midpoint,y_midpoint,"r-.",label = "Mid-point Approximation")
    plt.plot(x_heun,y_heun,"g--",label = "Heun Approximation")
    plt.plot(x_ralston,y_ralston,"k--",label = "Ralston Approximation")
    plt.xlabel("Time")
    plt.ylabel("Temperature")
    plt.legend()
    plt.show()
