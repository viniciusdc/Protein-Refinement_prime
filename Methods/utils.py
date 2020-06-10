import numpy as np
import random
from scipy.linalg import svd


# funções auxiliares utilizadas:
# prod_interno :: produto interno entre pares (A,v) com A matriz (n, p) e v vetor m-dimensional;
# def :: <(A, v), (B, w)> = Tr(AB') + sum_{i = 1}^{m} (v_i * w_i)
def prod_interno(matrix_a, matrix_b, v, w):
    # n :: número de linhas de A;
    # p :: número de linhas de B;
    n, p = matrix_a.shape
    if (n, p) != matrix_b.shape:
        print(' Dimensões incompativeis entre as matrizes !')
        return
    # m :: tamanho/dimensão de v;
    m = len(v)
    if m != len(w):
        print(' Dimensões incompativeis entre os vetores !')
        return

    prod_soma = 0.0
    for j in range(p):
        for i in range(n):
            prod_soma += matrix_a[i, j] * matrix_b[i, j]
    for i in range(m):
        prod_soma += v[i] * w[i]

    return prod_soma


# distância Euclidiana entre as linhas i e j da matriz X;
def distance(i, j, matrix_x):
    dist = 0.0
    # dim :: dimensão da matriz X;
    dim = matrix_x.shape[1]

    for k in range(dim):
        dist += (matrix_x[i, k] - matrix_x[j, k]) ** 2

    return np.sqrt(dist)


# centralizar um conjunto de pontos x;
def centralizar(x):
    # caso a entrada não esteja no formato adequado;
    x = np.array(x, dtype=float)

    # k :: dimensão
    # n :: número de pontos
    n, k = x.shape

    ponto_medio = (1 / n) * np.dot(np.ones(n), x)
    # ponto transladado;
    x = x - (np.ones((n, 1)) * ponto_medio)

    return x


# rmsd :: calcula a raiz quadrada da média dos desvios entre duas estruturas A e B (centralizadas);
def rmsd(matrix_a, matrix_b):
    # n :: número de pontos;
    n, k = matrix_a.shape
    if matrix_a.shape != matrix_b.shape:
        print('Dimensões incompativeis entre as matrizes')
        print('Dimensões da matriz A : {}'.format(matrix_a.shape))
        print('Dimensões da matriz B : {}'.format(matrix_b.shape))
        return 'NaN'
    # Procrustes:
    # Given two matrices A and B it is asked to find an orthogonal matrix Q which most closely maps A to B.
    matrix_a = centralizar(matrix_a)
    matrix_b = centralizar(matrix_b)
    singular_value_dec = svd(np.dot(matrix_b, matrix_a.T))
    matri_q = np.dot(singular_value_dec[0], singular_value_dec[2])

    correlation = np.linalg.norm(np.dot(matri_q, matrix_a) - matrix_b, ord='fro')

    return np.sqrt(1 / n) * correlation


# mde :: erro médio dos desvios;
# lb :: limitantes inferiores;
# ub :: limitantes superiores;
def mde(ponto, vec_u, vec_v, lb, ub):
    # dim :: tamanho/ dimensão do vetor u;
    # m :: número de limitantes/distâncias conhecidos/das;
    dim = len(vec_u)
    m = len(lb)

    soma = 0
    for s in range(dim):
        dist_ponto = distance(vec_u[s], vec_v[s], ponto)
        soma += max((lb[s] - dist_ponto) / lb[s], 0) + max((dist_ponto - ub[s]) / ub[s], 0)

    return soma / m


# apenas uma checagem para verificar se o problema de dimensionamento foi resolvido.
def check_solution_dimension(matrix_a, matrix_b):
    if matrix_a.shape[0] != matrix_b.shape[0]:
        return True
    else:
        return False


# dado um conjunto de pontos, projeta-se suas distâncias sobre os limitantes lb e ub;
# saída :: vetor armazenando as distâncias projetadas.
def dist_matrix_projection(dim, vec_u, vec_v, lower_bound, upper_bound, point_set, multistart=False):
    y = np.zeros(dim)
    if multistart:
        for s in range(dim):
            # d = distance(vec_u[s], vec_v[s], point_set)
            # if lower_bound[s] <= d <= upper_bound[s]:
            #     y[s] = d
            # elif d < lower_bound[s]:
            #     y[s] = random.uniform(lower_bound[s], upper_bound[s])
            # elif d > upper_bound[s]:
            #     y[s] = random.uniform(lower_bound[s], upper_bound[s])
            y[s] = random.uniform(lower_bound[s], upper_bound[s])
    else:
        for s in range(dim):
            d = distance(vec_u[s], vec_v[s], point_set)
            if lower_bound[s] <= d <= upper_bound[s]:
                y[s] = d
            elif d < lower_bound[s]:
                y[s] = lower_bound[s]
            elif d > upper_bound[s]:
                y[s] = upper_bound[s]
    return y
