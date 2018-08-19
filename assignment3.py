'''=======================================
! Lab No: 02
! Title : LU Decomposition
! Date: 13/08/2018
! Name : Ashutosh Kumar Mandal
! Roll No: 160121010
! Webmail : ashutosh3004@iitg.ac.in
=========================================='''


import numpy as np

#Creates the pivoting matrix for A. 
def pivotize(A):
    n = len(A)
    ID = np.identity(n)
    for k in range(n-1):
        maxindex = abs(A[k:,k]).argmax() + k
        if A[maxindex, k] == 0:
            raise ValueError("Matrix is singular.")
        # Swap
        if maxindex != k:
            ID[[k,maxindex]] = ID[[maxindex, k]]
    return ID

#Returns L, U and Pivot matrix for matrix A.
def lu(A):
    n = len(A)
    L = np.identity(n)
    U = np.zeros((n,n))
    P = pivotize(A)
    A2 = np.matmul(P, A)
    
    # Elimination
    for k in range(n-1):
        for row in range(k+1, n):
            multiplier = A2[row,k]/A2[k,k]
            A2[row, k:] = A2[row, k:] - multiplier*A2[k, k:]
            L[row,k] = multiplier
    U = A2
    return (L, U, P)

#Forward substitution function for lower triangular matrix
def forward(A, b):
    n = len(A)
    x = np.zeros(n)
    sum = 0
    for k in range(0, n):
        for j in range(k):
            sum = sum + x[j]*A[k][j]
        x[k] = (b[k] - sum)/A[k][k]
       
    return x

#Backward substitution function for upper triangular matrix
def backward(U, b):
    n = len(U)
    x = np.zeros(n)
    x[n-1] = b[n-1]/U[n-1][n-1]
    for k in range(n-2, -1, -1):
        temp = b[k]
        for j in range(k+1, n):
            temp = temp - x[j]*U[k][j]
        x[k] = temp/U[k][k]
        
    return x

#Function to create inverse for upper triangular matrix
def uinv(U):
    n = len(U)
    I = np.identity(n)
    #Eliminate Upper part
    for k in range(n-1, -1, -1):
        for row in range(k-1, -1, -1):
            multiplier = U[row,k]/U[k,k]
            U[row, k:] = U[row, k:] - multiplier*U[k, k:]
            I[row, k:] = I[row, k:] - multiplier*I[k, k:]
    
    #Normalize diagonal elements
    for k in range(n):
        I[k] = I[k]/U[k][k]
        U[k] = U[k]/U[k][k]
    return (I, U)

#Function to create inverse for lower triangular matrix L whose
            #diagonal elements are 1.
def linv(L):
    n = len(L)
    for k in range(n):
        for j in range(n):
            if (k != j):
                L[k][j] = -1 * L[k][j]
    return L

#Function to solve Ax = b using LU decomposition
#Returns the solution
def lusolve(A, b):
    n = len(A)
    (L, U, P) = lu(A)
    Y = forward(L, b)
    X = backward(U, Y)
    return X
           

#Driver code
if __name__ == "__main__":
    print(" Enter input csv file (e.g. file.csv) : \n") 
    file = input()
    #Read input from csv file
    data = np.genfromtxt(file, delimiter=',')
    
    #Separate values into matrix A and vector b
    A = data[:,0:-1]
    print("\n\n Matrix A : \n ",A)
    b = data[:,-1]
    
    #print("\n\n Solutions : \n",lusolve(A, b))
    
    (L, U, P) = lu(A)
    print("\n\n Pivot matrix : \n", P, "\n\n L : \n", L, "\n\n U : \n", U)
    print("\n\n Product of L and U : \n",np.matmul(L, U))
    
    y = forward(L, b)
    print("\n\n Y = UX = : \n", forward(L, b))
    x = backward(U, y)
    print("\n\n Solution = X : \n", x)
    
    (invU, modU) = uinv(U)
    invL = linv(L)
    print("\n\n Inverse of L : \n", invL, "\n\n Inverse of U : \n", invU)
    print("\n\n Inverse of A : \n", np.matmul(invU,invL)) #The inverse calculated will be wrong for a matrix requiring pivoting.
    