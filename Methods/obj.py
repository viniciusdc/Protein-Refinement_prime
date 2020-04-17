import numpy as np
from Methods.utils import distance


# stress
def stress(x, y, u, v, w):
    m = len(u)
    soma = 0.0

    for k in range(m):
        prod = distance(u[k], v[k], x) - y[k]
        soma = soma + w[k] * (prod * prod)
    return soma


# stress gradient
def grad_stress(x, y, u, v, w):
    x = np.array(x)

    (n, p) = x.shape
    m = len(y)
    gx = np.zeros(x.shape)
    gy = np.zeros(m)
    memory = np.zeros(n)

    # Calculation for Gx and Gy (compact form, except the case when u == v)
    for k in range(m):
        temp = distance(u[k], v[k], x)
        gy[k] = - 2.0 * w[k] * (temp - y[k])
        if temp > 0.0:
            temp = - w[k] * (y[k] / temp)
            memory[u[k]] = memory[u[k]] + w[k] + temp
            memory[v[k]] = memory[v[k]] + w[k] + temp
            temp = 2.0 * (-w[k] - temp)
            for j in range(p):
                gx[u[k], j] += temp * x[v[k], j]
                gx[v[k], j] += temp * x[u[k], j]

    # case when u == v
    for j in range(p):
        for i in range(n):
            gx[i, j] += 2.0 * memory[i] * x[i, j]

    return gx, gy
