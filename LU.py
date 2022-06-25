import scipy
import scipy.linalg # SciPy Linear Algebra Library to compare time of executions ---> our algorithm vs scipy
import numpy as np
import pandas as pd #to upload the data on dataFrame and save it as csv
import time #to measure time of executions
import sys 

def validate(A):
     
     if len(A)!=len(A[0]):
        print("matrix is not square")
        return False
     for i in range(len(A)):
        if A[i,i]==0:
            print("matrix has zero in the main diagonal")
            return False

def LUdecomposition_Optimized(A):
    if validate(A)==False:
        sys.exit(0)
    #initialize n with the number of rows / columns
    n = A.shape[0]
    #initialize matrix U size: ---> nxn with all zeros ---> the zero matrix
    U = np.zeros((n, n), dtype=np.double)
    # initialize matrix L size: ---> nxn to be the identity matrix ---> ones in the main diagonal
    L = np.eye(n, dtype=np.double)
    
    for k in range(n):
        
        U[k, k:] = A[k, k:] - np.matmul(L[k,:k] , U[:k,k:])
        L[(k+1):,k] = (A[(k+1):,k] - np.matmul(L[(k+1):,:] , U[:,k])) / U[k, k]
       
    
    return L, U






def luDecomposition_NonOptimized(A):
    
    if validate(A)==False:
       sys.exit(0)
    #initialize n with the number of rows / columns
    n = A.shape[0]
    
    #initialize matrices lower ,upper size nxn

    lower = [[0 for x in range(n)]
                 for y in range(n)]
    upper = [[0 for x in range(n)]
                 for y in range(n)]
    
   # Decomposing matrix into Upper
    # and Lower triangular matrix
    for i in range(n):
 
        # Upper Triangular
        for k in range(i, n):
 
            # Summation of L(i, j) * U(j, k)
            sum = 0.0
            for j in range(i):
                sum += (lower[i][j] * upper[j][k])
 
            # Evaluating U(i, k)
            upper[i][k] = A[i][k] - sum
 
        # Lower Triangular
        for k in range(i, n):
            if (i == k):
                lower[i][i] = 1  # Diagonal as 1
            else:
 
                # Summation of L(k, j) * U(j, i)
                sum = 0.0
                for j in range(i):
                    sum += (lower[k][j] * upper[j][i])
 
                # Evaluating L(k, i)
                lower[k][i] = float((A[k][i] - sum) /
                                  upper[i][i])
    return lower,upper





 



    """
 ____________  main-test _____________________________
"""

print("start reading file...") 

fileName = "data_10.txt"

input_matrix = np.loadtxt(fileName, delimiter=',' , dtype = np.float64)



print("data fetched successfully...")
print("number of lines in file:", input_matrix.shape[0])
while(True):
     choice = int(input("Enter which version do you want to use? \n 0 - non-optimized \n 1 - optimized\n"))
     if choice == 0 or choice == 1:
         break 
        
     else:
         print("\n invalid option please enter again...")

print("LU decomposition start...")
#measuring time using the python library
print("starting scipy python LU decomposition...")
t0 = time.time()
scipy.linalg.lu(input_matrix) # --> python scipy LU decomposition
t1 = time.time() - t0


print("python scipy LU decomposition time:",t1)
#measuring time using our LU decomposition code
print("starting OUR LU decomposition...")
t2 = time.time()

if choice == 1:
   lower , upper = LUdecomposition_Optimized(input_matrix) 
elif choice == 0:
       lower , upper = luDecomposition_NonOptimized(input_matrix) 
t3 = time.time() - t2

print("Our LU decomposition time:",t3)
#to check if input matrix is the same as LU
res = np.matmul(lower,upper);
for i in range(res.shape[0]):
  for j in range(res.shape[0]):
        if (res[i][j] != input_matrix[i][j]):
         isIdentical=0
        else:
         isIdentical=1

if(isIdentical == 1):
    print("matrix identical")
else:
    print("matrix not identical")
    
#to check if input matrix is the same as LU
print(np.allclose(input_matrix, np.matmul(lower,upper)))

print("saving to file...") 
t4 = time.time()
print("saving to file lower...") 
DFL = pd.DataFrame(lower)
DFL.to_csv("L.csv")
print("saving to file upper...") 

DFU = pd.DataFrame(upper)
DFU.to_csv("U.csv")

arrL = np.array(lower)
arrU = np.array(upper)

np.savetxt('L.txt',arrL, delimiter=',',fmt="%.4f")
np.savetxt('U.txt',arrU, delimiter=',',fmt="%.4f")

t5 = time.time() - t4
print("time to save to files: ", t5)
