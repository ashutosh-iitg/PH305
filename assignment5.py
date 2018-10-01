'''======================================
! Lab No: 05
! Title : Power Method
! Date: 10/09/2018
! Name : Ashutosh Kumar Mandal
! Roll No: 160121010
! Webmail : ashutosh3004@iitg.ac.in
========================================='''

import numpy as np
from matplotlib import pyplot as plt


if __name__ == '__main__':
    np.set_printoptions(formatter={'float': lambda x: "{0:0.6f}".format(x)})
    A = np.genfromtxt('file05.csv', delimiter = ',')
    x_pre = [1, 2, 3]
    x_curr = [1, 2, 3]
    print(" Matrix A : \n ", A,"\n\n X_Previous : \n ",x_pre,"\n\n X_Current : \n ",x_curr, "\n\n")
    error = 1
    n_iter = 0
    
    f =  open("out1.txt",'a')
    while error > 0.0000000001:
        n_iter += 1
        x_pre[:] = x_curr[:]
        x_curr = np.matmul(A,x_pre)
        max_val = x_curr[np.argmax(x_curr)]
        x_curr = x_curr/max_val
        error = np.linalg.norm(x_curr - x_pre)
        print(n_iter, " ", "{:.6f}".format(max_val), " ", "{:.6f}".format(error), file = f)
        #f.write(str(n_iter) + " " + str(max_val.round(6)) + " " + str(error.round(6)) + "\n")
        print(" ", "{:.6f}".format(error))
    print("\n Largest Eigen Value : ", "{:.6f}".format(max_val))
    print(" Corresponding Eigen Vector : ", x_curr)
    f.close()
    
    inverse = np.linalg.inv(A)
    x_pre = [1, 2, 3]
    x_curr = [1, 2, 3]
    print("\n\n\n")
    error = 1
    n_iter = 0
    
    g = open("out2.txt",'a')
    while error > 0.0000000001:
        n_iter += 1
        x_pre[:] = x_curr[:]
        x_curr = np.matmul(inverse,x_pre)
        max_val = x_curr[np.argmax(x_curr)]
        x_curr = x_curr/max_val
        error = np.linalg.norm(x_curr - x_pre)
        print(n_iter, " ", "{:.6f}".format(1/max_val), " ", "{:.6f}".format(error), file = g)
        #g.write(str(n_iter) + " " + str((1/max_val).round(6)) + " " + str(error.round(6)) + "\n")
        print(" ", "{:.6f}".format(error))
    print("\n Smallest Eigen Value : ","{:.6f}".format(1/max_val))
    print(" Corresponding Eigen Vector : ", x_curr)
    g.close()
    
    #Define matrix A
    A = np.zeros((10,10))
    for i in range(10):
        for j in range(10):
            if i == j:
                A[i][i] = 4
            elif (abs(i-j) == 1):
                A[i][j] = 2
            elif (abs(i-j) == 2):
                A[i][j] = 1
    print(A)
    
    x_pre = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    x_curr = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    print(" Matrix A : \n ", A,"\n\n X_Previous : \n ",x_pre,"\n\n X_Current : \n ",x_curr, "\n\n")
    error = 1
    n_iter = 0
    
    f =  open("out3.txt",'a')
    while error > 0.0000000001:
        n_iter += 1
        x_pre[:] = x_curr[:]
        x_curr = np.matmul(A,x_pre)
        max_val = x_curr[np.argmax(x_curr)]
        x_curr = x_curr/max_val
        error = np.linalg.norm(x_curr - x_pre)
        print(n_iter, " ", "{:.6f}".format(max_val), " ", "{:.6f}".format(error), file = f)
        #f.write(str(n_iter) + " " + str(max_val.round(6)) + " " + str(error.round(6)) + "\n")
        print(" ", "{:.6f}".format(error))
    print("\n Largest Eigen Value : ", "{:.6f}".format(max_val))
    print(" Corresponding Eigen Vector : ", x_curr)
    f.close()
    
    inverse = np.linalg.inv(A)
    x_pre = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    x_curr = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    print("\n\n\n")
    error = 1
    n_iter = 0
    
    g = open("out4.txt",'a')
    while error > 0.0000000001:
        n_iter += 1
        x_pre[:] = x_curr[:]
        x_curr = np.matmul(inverse,x_pre)
        max_val = x_curr[np.argmax(x_curr)]
        x_curr = x_curr/max_val
        error = np.linalg.norm(x_curr - x_pre)
        print(n_iter, " ", "{:.6f}".format(1/max_val), " ", "{:.6f}".format(error), file = g)
        #g.write(str(n_iter) + " " + str((1/max_val).round(6)) + " " + str(error.round(6)) + "\n")
        print(" ", "{:.6f}".format(error))
    print("\n Smallest Eigen Value : ","{:.6f}".format(1/max_val))
    print(" Corresponding Eigen Vector : ", x_curr)
    g.close()
    print(np.linalg.eig(A))
