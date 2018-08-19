'''=======================================
! Lab No: 02
! Title : Gauss Jordan Elimination
! Date: 13/08/2018
! Name : Ashutosh Kumar Mandal
! Roll No: 160121010
! Webmail : ashutosh3004@iitg.ac.in
=========================================='''

import numpy as np


def GaussJ(A, B, b, doPricing = True):
    '''
    Gauss Jordan elimination with partial pivoting.
    
    input: A is an n x n numpy matrix
           B is an n x n identity matrix
           b is an n x 1 numpy array
    output: x is the solution of Ax=b 
            with the entries permuted in 
            accordance with the pivoting 
            done by the algorithm
    post-condition: A, B and b have been modified.
                    B is modified to be the inverse of A.
    '''
    n = len(A)
    if b.size != n:
        raise ValueError("Invalid argument: incompatible sizes between"+
                         "A & b.", b.size, n)
    # k represents the current pivot row. Since GE traverses the matrix in the 
    # upper right triangle, we also use k for indicating the k-th diagonal 
    # column index.
    
    # Elimination
    for k in range(n-1):
        if doPricing:
            # Pivot
            maxindex = abs(A[k:,k]).argmax() + k
            if A[maxindex, k] == 0:
                raise ValueError("Matrix is singular.")
            # Swap
            if maxindex != k:
                A[[k,maxindex]] = A[[maxindex, k]]
                B[[k,maxindex]] = B[[maxindex, k]]
                b[[k,maxindex]] = b[[maxindex, k]]
        else:
            if A[k, k] == 0:
                raise ValueError("Pivot element is zero. Try setting doPricing to True.")
        #Eliminate lower part
        for row in range(k+1, n):
            multiplier = A[row,k]/A[k,k]
            A[row, k:] = A[row, k:] - multiplier*A[k, k:]
            B[row, k:] = B[row, k:] - multiplier*B[k, k:]
            b[row] = b[row] - multiplier*b[k]
    
    #Eliminate Upper part
    for k in range(n-1, -1, -1):
        for row in range(k-1, -1, -1):
            multiplier = A[row,k]/A[k,k]
            A[row, k:] = A[row, k:] - multiplier*A[k, k:]
            B[row, k:] = B[row, k:] - multiplier*B[k, k:]
            b[row] = b[row] - multiplier*b[k]
    
    #Normalize diagonal elements
    for k in range(n):
        B[k] = B[k]/A[k][k]
        b[k] = b[k]/A[k][k]  
        A[k] = A[k]/A[k][k]
            
    x = np.zeros(n)
    x = b
    
    return x



if __name__ == "__main__":
    
    print(" Enter input csv file (e.g. file.csv) : \n") 
    file = input()
    
    #Read input from csv file
    data = np.genfromtxt(file, delimiter=',')
    
    #Separate values into matrix A and vector b
    A = data[:,0:-1]
    b = data[:,-1]
    
    #Create an identity matrix of size of A
    B = np.identity(len(A))
    
    #print(GaussJ(np.copy(A), np.copy(B), np.copy(b), doPricing = False))
    print("Solutions : \n", GaussJ(A,B,b))
    #print(A)
    print("\n Inverse of A : \n", B)
