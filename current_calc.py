

BOLTZMANN_CONSTANT_J = 1.380649 * pow(10,-23)

import numpy as np

def beta(temperature):
    beta = 1/(temperature * BOLTZMANN_CONSTANT_J)
    return beta


def A(x, y, n):
    A = (x * (7 + 2 * (x ** 2) * (n ** 2)) * (1 - np.cos(y * n))) / (
                (1 + (x ** 2) * (n ** 2)) * (16 + (x ** 2) * (n ** 2)))
    return A


def B(x, y, n, m):
    num1 = x * (5 * n + (n + m) * (4 - (x ** 2) * (n ** 2)))
    num2 = np.sin(y * n) * (1 - np.cos(y * m)) + np.sin(y * m) * (1 - np.cos(y * n))
    denom = (n + m) * m * (1 + (x ** 2) * (n ** 2)) * (16 + (x ** 2) * (n ** 2)) * (1 + (x ** 2) * ((n + m) ** 2))
    B = num1 * num2 / denom
    return B

def C(x, y):
    C = 0
    n = 20
    for i in range(n):
        C += A(x, y, i)
    return C

def D(x, y):
    D = 0
    n = 20
    m = 20
    for i in range(n):
        i_for_j = i + 1
        i_fixed = i + 2
        for j in range(i_for_j):
            j_fixed = -1 * (j + 1)
            D += (B(x, y, i_fixed, j_fixed) + B(x, y, j_fixed, i_fixed))

    for k in range(n):
        k_fixed = k + 1
        D += B(x, y, k_fixed, 1)

    return D

def G(x,y, delta):
    G = (1 - (2 * delta)) * C(x, y) - ((1 / (2 * np.pi))) * D(x, y)
    return G

def phi_1(xi, delta, x, y):
    phi_1 = (12 / (np.pi ** 2)) * x * C(x, y)
    return phi_1

def phi_2(xi, delta, x, y):
    phi_2 = (12 / (np.pi ** 2)) * x * G(x, y, delta)
    return phi_2

def get_velocity(period, L, diffusion, a1, a2, alpha, temperature, dc):
    xi = ((L/(2 * np.pi)) ** 2) / (diffusion * period)
    x = (2 * np.pi) * xi
    # tau_ = 0.6
    # tau = 1
    delta = dc
    y = (2 * np.pi) * delta
    num1 = diffusion * np.pi * (beta(temperature) ** 3) * (a1 ** 2) * a2 * ((1 - alpha) ** 2)
    num2 = ((1 + alpha) * phi_1(xi, delta, x, y)) + ((1 - alpha) * phi_2(xi, delta, x, y))
    denom = (4 * L)
    velocity = (num1 * num2) / (denom)
    return velocity

def get_current(velocity, ne, sigma, q):
    current = velocity * ne * sigma * q
    return current