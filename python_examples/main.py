
import numpy.matlib as np
import scipy as sp
import time

N = 128
# A = [[1, 2, 3], [0, 2, 5], [0, 0, 9]]

def gen_rand(N):
    A = [[0] * N for i in range(N)]
    for i in range(N):
        for j in range(N):
            A[i][j] = np.random.rand()
    return A

A = gen_rand(N)
start_time1 = time.time()
p, l, u = sp.linalg.lu(A)
print(time.time() - start_time1)
# print(matrix)
# print(l * u)

def matrix_mult(A, B, N):
    C = [[0] * N for i in range(N)]
    for i in range(N):
        for j in range(N):
            C[i][j] = 0
            for k in range(N):
                C[i][j] += A[i][k] * B[k][j]
    # print(C)

def lu(A, N):
    L = [[0] * N for i in range(N)]
    for i in range(N):
        L[i][i] = 1
        for j in range(i + 1, N):
            L[j][i] = A[j][i] / A[i][i]
            A[j][i] = 0
            for k in range(i, N):
                A[j][k] -= A[i][k] * L[j][i]
    # print(L)
    # print(A)
    matrix_mult(L, A, N)

def lu_np(A, N):
    L = np.zeros((N, N))
    for i in range(N):
        L[i][i] = 1
        # L[i+1:][i] = (1 / A[i][i]) * A[i+1:][i]
        for j in range(i + 1, N):
            L[j][i] = A[j][i] / A[i][i]
            A[j][i] = 0
            A[j][i+1:] -= L[j][i] * A[i][i+1:]
    # print(L)
    # print(A)
    matrix_mult(L, A, N)

start_time2 = time.time()
lu(A, N)
print(time.time() - start_time2)

A = gen_rand(N)
start_time3 = time.time()
lu_np(A, N)
print(time.time() - start_time3)
