'''=======================================
! Lab No: 03
! Title : Cholesky Decomposition
! Date: 20/08/2018
! Name : Ashutosh Kumar Mandal
! Roll No: 160121010
! Webmail : ashutosh3004@iitg.ac.in
==========================================='''

import numpy as np

def cholesky(A):
    """Decomposes a symmetric matrix A such that A = LL'.
       Input: A symmetric matrix A
       Output: A lower triangular matrix L"""
    n = len(A)
    L = np.zeros((n,n))
    for i in range(n):
        for j in range(i+1):
            sum = 0;
            if (j == i):
                for k in range(j):
                    sum += L[j][k] ** 2
                L[j][j] = np.sqrt(A[j][j] - sum)
            else:
                for k in range(j):
                    sum += L[i][k]*L[j][k]
                L[i][j] = (A[i][j] - sum) / L[j][j]       
    return L;

#Function to create transpose of a square matrix
def transpose(M):
    n = len(M)
    trans_M = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            trans_M[j][i] = M[i][j]
    return trans_M
    

if __name__ == "__main__":
    
    print(" Enter input csv file (e.g. file.csv) : \n") 
    file = input()
    
    #Read input from csv file
    A = np.genfromtxt(file, delimiter=',')
    L = cholesky(A)
    trL = transpose(L)
    print("\n\n A : \n", A)
    print("\n\n L : \n", L)
    print("\n\n L Transpose : \n",trL)
    print("\n\n Product of L and L transpose : \n", np.matmul(L,trL))