import numpy as np

N = 64
np.random.seed(0)

# Constant diagonals
dp  = np.complex64(5 + 0.5j)
dp1 = np.complex64(5.5 + 0.5j)
dp2 = np.complex64(4.5 + 0.5j)
do  = np.complex64(1 - 0.2j)

# Right hand side vector
b = (np.random.randn(N) + 1j*np.random.randn(N)).astype(np.complex64)

# Build matrix and solve for reference
A = np.zeros((N, N), dtype=np.complex64)
A[0,0] = dp1
A[-1,-1] = dp2
for i in range(1, N-1):
    A[i,i] = dp
for i in range(N-1):
    A[i+1,i] = do
    A[i,i+1] = do

x = np.linalg.solve(A, b)

def save_vector(fname, arr):
    with open(fname, 'w') as f:
        for v in arr:
            f.write(f"{v.real} {v.imag}\n")

save_vector('b.dat', b)
save_vector('golden.dat', x)

print('data generated')
