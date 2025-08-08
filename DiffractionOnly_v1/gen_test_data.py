import numpy as np
from pathlib import Path

N = 64
Lx = 45e-6
Ly = 45e-6
dx = Lx / N
dy = Ly / N
dz = 1e-4
k = 7853981.6339
eps = 1e-12


def custom_thomas_solver(dp, dp1, dp2, do, b):
    n = b.shape[0]
    c_prime = np.zeros(n-1, dtype=np.complex64)
    d_prime = np.zeros(n, dtype=np.complex64)
    c_prime[0] = do / dp1
    d_prime[0] = b[0] / dp1
    for i in range(1, n-1):
        denom = dp - do * c_prime[i-1]
        c_prime[i] = do / denom
        d_prime[i] = (b[i] - do * d_prime[i-1]) / denom
    d_prime[n-1] = (b[n-1] - do * d_prime[n-2]) / (dp2 - do * c_prime[n-2])
    x = np.zeros(n, dtype=np.complex64)
    x[n-1] = d_prime[n-1]
    for i in range(n-2, -1, -1):
        x[i] = d_prime[i] - c_prime[i] * x[i+1]
    return x


def compute_b_vector(dp, dp1, dp2, do, x0):
    n = x0.shape[0]
    b = np.zeros(n, dtype=np.complex64)
    b[0] = dp1 * x0[0] + do * x0[1]
    for i in range(1, n-1):
        b[i] = do * x0[i-1] + dp * x0[i] + do * x0[i+1]
    b[n-1] = do * x0[n-2] + dp2 * x0[n-1]
    return b


def adi_x(phi):
    ung = np.complex64(1j * dz / (4 * k * dx * dx))
    phi_inter = np.zeros_like(phi, dtype=np.complex64)
    for j in range(N):
        if abs(phi[1, j]) < eps:
            ratio_x0 = np.float32(1.0)
        else:
            ratio_x0 = phi[0, j] / phi[1, j]
        if abs(phi[N-2, j]) < eps:
            ratio_xn = np.float32(1.0)
        else:
            ratio_xn = phi[N-1, j] / phi[N-2, j]
        dp_B = np.float32(1.0) - 2 * ung
        dp1_B = dp_B + ung * ratio_x0
        dp2_B = dp_B + ung * ratio_xn
        b = compute_b_vector(dp_B, dp1_B, dp2_B, ung, phi[:, j])
        dp_A = np.float32(1.0) + 2 * ung
        dp1_A = dp_A - ung * ratio_x0
        dp2_A = dp_A - ung * ratio_xn
        phi_inter[:, j] = custom_thomas_solver(dp_A, dp1_A, dp2_A, -ung, b)
    return phi_inter


def adi_y(phi):
    ung = np.complex64(1j * dz / (4 * k * dy * dy))
    phi_inter = np.zeros_like(phi, dtype=np.complex64)
    for i in range(N):
        if abs(phi[i, 1]) < eps:
            ratio_y0 = np.float32(1.0)
        else:
            ratio_y0 = phi[i, 0] / phi[i, 1]
        if abs(phi[i, N-2]) < eps:
            ratio_yn = np.float32(1.0)
        else:
            ratio_yn = phi[i, N-1] / phi[i, N-2]
        dp_B = np.float32(1.0) - 2 * ung
        dp1_B = dp_B + ung * ratio_y0
        dp2_B = dp_B + ung * ratio_yn
        b = compute_b_vector(dp_B, dp1_B, dp2_B, ung, phi[i, :])
        dp_A = np.float32(1.0) + 2 * ung
        dp1_A = dp_A - ung * ratio_y0
        dp2_A = dp_A - ung * ratio_yn
        phi_inter[i, :] = custom_thomas_solver(dp_A, dp1_A, dp2_A, -ung, b)
    return phi_inter


def save_matrix(fname, arr):
    with open(fname, 'w') as f:
        for j in range(N):
            for i in range(N):
                v = arr[i, j]
                f.write(f"{float(v.real)} {float(v.imag)}\n")


def main():
    np.random.seed(0)
    phi_in = (np.random.randn(N, N) + 1j*np.random.randn(N, N)).astype(np.complex64)
    tmp = adi_x(phi_in)
    phi_out = adi_y(tmp)
    here = Path(__file__).resolve().parent
    save_matrix(here / 'phi_in.dat', phi_in)
    save_matrix(here / 'golden.dat', phi_out)
    print('data generated')

if __name__ == '__main__':
    main()
