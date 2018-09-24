import numpy as np
import math
import copy

np.set_printoptions(formatter={'float': lambda x: "{0:0.6f}".format(x)})

def max_elem(A,n):
    max_val = 0
    p = 0
    q = 0
    for i in range(n):
        for j in range(i):
            if (abs(max_val) < abs(A[i][j])):
                max_val = A[i][j]
                p = j
                q = i
    return max_val,p,q

def symm(A,n):
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            elif A[i][j] != A[j][i]:
                raise ValueError("Matrix is not symmetric.")
    return 1
    
def rotate(A,n,p,q):
    R = np.identity(n)
    D = copy.deepcopy(A)
    theta = (A[q][q] - A[p][p])/(2 * A[p][q])
    t = np.sign(theta)/(abs(theta) + math.sqrt(theta**2 + 1))
    c = 1/math.sqrt(t**2 + 1)
    s = c * t
    
    R[p][p] = c
    R[p][q] = s
    R[q][p] = -s
    R[q][q] = c
    
    #D = np.matmul(np.matmul(np.transpose(R),A),R)
    
    D[p][q] = 0
    D[q][p] = 0
    D[p][p] = c**2 * A[p][p] + s**2 * A[q][q] - 2*c*s * A[p][q]
    D[q][q] = s**2 * A[p][p] + c**2 * A[q][q] + 2*c*s * A[p][q]
    for i in range(n):
        if ((i != p) and (i != q)):
            D[i][p] = c * A[i][p] - s * A[i][q]
            D[p][i] = D[i][p]
            D[i][q] = c * A[i][q] + s * A[i][p]
            D[q][i] = D[i][q]
    
    return D,R

def jacobi(A,n):
    D = copy.deepcopy(A)
    tol = 1e-10
    R_prev = np.identity(n)
    n_iter = 0
    
    for i in range(10*(n*n)):
        n_iter += 1
        (m,p,q) = max_elem(A,n)
        if (abs(m) < tol):
            break
        (A,R_curr) = rotate(A,n,p,q)
        R_prev = np.matmul(R_prev,R_curr)
        
    print("\n\n Diagonalized Matrix :\n",A)
    print("\n\n Transformation Matrix :\n",R_prev)
    print("\n\n Number of iterations taken : ",n_iter)
    return A
    

if __name__ == "__main__" :
    
    A = np.genfromtxt('file06.csv',delimiter = ',')
    print("\n Input Matrix :\n",A)
    n = len(A)
    #print(np.linalg.eig(A))
    
    if(symm(A,n) == 1):
        A = jacobi(A,n)