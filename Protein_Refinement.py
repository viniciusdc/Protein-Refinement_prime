import numpy as np
from scipy.linalg import svd
import time
from Methods.distance_file_gen import gen_distance_file

# recebe como input o nome da proteína a ser refinada

print(">> Digite o nome da proteína: ")
filename = str(input())

# Arquivos guia:

' -- Gerador do arquivo de distâncias -- '
# raid :: root archive index directory
local_dir = 'C:\\Users\\viniv\\Desktop\\Testes'

# - Arquivo PDB:
# raid_pdb :: diretório do arquivo PDB principal.
raid_pdb = local_dir + '\\Teste {}\\{}.txt'.format(filename, filename)
pdb = np.genfromtxt(raid_pdb, dtype='str')

print(":: Arquivo PDB lido com sucesso ! Iniciando gerador de distâncias;")

gen_distance_file(pdb, filename, local_dir)

print(":: Processo finalizado com êxito, aguardando leitura dos dados...")

# abertura do arquivo de distâncias gerado;
distancias = np.genfromtxt(raid_d, dtype='str')

# variáveis de ajuste:
# Noise :: Grau de perturbação da solução esperada;
# TOL :: Tolerância para o SPG;
# N :: Número máximo de iterações aceito;
# M :: parâmetro de não monotonia de GLL;
Noise = 1e0
TOL = 1e-10
N = 2000
M = 15

# Informações sobre a proteína:

" -- Funções auxiliares -- "


def prod_interno(A1, B1, A2, B2):
    n, p = A1.shape
    m = len(A2)
    Prod = 0.0
    for j in range(p):
        for i in range(n):
            Prod += A1[i, j] * B1[i, j]
    for i in range(m):
        Prod += A2[i] * B2[i]

    return Prod


def distance(a, b, X):
    dist = 0.0
    p = len(X[0])
    for i in range(p):
        dist += (X[a, i] - X[b, i]) ** 2
    return np.sqrt(dist)


def centralizar(x):
    x = np.array(x, dtype=float)
    Mean = np.array(x[0, :])
    (m, n) = x.shape

    for i in range(1, m):
        Mean += x[i, :]
    Mean = Mean / m
    for i in range(m):
        x[i, :] += - Mean

    return x


" -- Função objetivo e gradiente -- "


def stress(x, y, u, v, w):
    m = len(u)
    soma = 0.0

    for k in range(m):
        prod = distance(u[k], v[k], x) - y[k]
        soma = soma + w[k] * (prod * prod)
    return soma


def grad_stress(x, y, u, v, w):
    x = np.array(x)

    (n, p) = x.shape
    m = len(y)
    Gx = np.zeros(x.shape)
    Gy = np.zeros(m)
    memory = np.zeros(n)

    # Cálculo de Gx e Gy (forma compacta, todos os passos em um, exceto o caso u == v)
    for k in range(m):
        temp = distance(u[k], v[k], x)
        Gy[k] = - 2.0 * w[k] * (temp - y[k])
        if temp > 0.0:
            temp = - w[k] * (y[k] / temp)
            memory[u[k]] = memory[u[k]] + w[k] + temp
            memory[v[k]] = memory[v[k]] + w[k] + temp
            temp = 2.0 * (-w[k] - temp)
            for j in range(p):
                Gx[u[k], j] += temp * x[v[k], j]
                Gx[v[k], j] += temp * x[u[k], j]

    # Completando o cálculo de Gx (caso u == v)
    for j in range(p):
        for i in range(n):
            Gx[i, j] += 2.0 * memory[i] * x[i, j]

    return Gx, Gy


# Método do Gradiente Projetado Espectral:


def spg(f, grad_f, Xo, yo, u, v, w, lb, ub, TOL, N, M):
    X = Xo
    Xp = Xo
    Yp = yo
    Y = yo
    Fo = f(Xo, yo, u, v, w)
    print('Valor inicial da função objetivo: {}'.format(Fo))
    Fn = Fo
    F = np.zeros(M)
    F[0] = Fo
    Gx, Gy = grad_f(Xo, yo, u, v, w)
    Gxp = Gx
    Gyp = Gy
    BackTracking = 0
    GtD = 0.0

    alpha = 1.0
    pB = 1.0
    k = 0
    print('Iteração máxima: {}, TOL: {}, M: {} e Ruído: {}.'.format(N, TOL, M, Noise))

    while k < N:
        if Fo < TOL:
            print('Valor da função objetivo menor do que a TOL: {}.'.format(TOL))
            return Xo, BackTracking, k, Fo / 2, GtD

        if k == 0:
            pB = 1.0
        if k != 0:
            Yx = Gx - Gxp
            Yy = Gy - Gyp
            Zx = X - Xp
            Zy = Y - Yp
            pB = prod_interno(Yx, Zx, Yy, Zy) / prod_interno(Zx, Zx, Zy, Zy)
            if pB < 1e-16:
                pB = 1e-16
            if pB > 1e16:
                pB = 1e16

        sX = Xo - (Gx / pB)
        sY = yo - (Gy / pB)

        for i in range(len(yo)):
            if sY[i] < lb[i]:
                sY[i] = lb[i]
            elif sY[i] > ub[i]:
                sY[i] = ub[i]

        DX = sX - Xo
        DY = sY - yo
        GtD = prod_interno(Gx, DX, Gy, DY)

        if abs(GtD) < 1e-10:
            print('Produto gradiente direção menor que a Tolerância {}'.format(1e-10))
            # print(' k: {:<3}, Fo: {:<10}, Fn: {:<10}, GtD: {:<10}, alpha: {: < 10}, pB: {:<10}, '
            #       .format(k, Fo / 2, Fn / 2, GtD, alpha, pB))
            return Xo, BackTracking, k, Fo / 2, GtD

        alpha = 1.0
        Xp = Xo
        Yp = yo
        Gxp = Gx
        Gyp = Gy

        X = Xp + alpha * DX
        Y = Yp + alpha * DY
        Fn = f(X, Y, u, v, w)
        F_max = max(F)

        while Fn > F_max + 1e-4 * alpha * GtD:
            BackTracking += 1
            alpha = alpha / 2
            X = Xp + alpha * DX
            Y = Yp + alpha * DY
            Fn = f(X, Y, u, v, w)

        # print('It: {:<3} Obj: {:<26} GtD: {:<26} pB: {:<26} alpha: {:<10}'.format(k, Fo / 2, GtD, pB, alpha))

        F[k % M] = Fn
        Xo = X
        yo = Y
        Fo = Fn
        Gx, Gy = grad_f(X, Y, u, v, w)
        k = k + 1

    return Xo, BackTracking, k, Fo / 2, GtD


# Pré-Inicialização:
" -- Reordenação dos átomos de acordo com o arquivo de distâncias -- "

# tomada dos átomos para ponto inicial Xo
L = max(max(np.array(distancias[:, 0], dtype='int')), max(np.array(distancias[:, 1], dtype='int')))

# Número de átomos da proteína, considerados após reordenação e exclusão.
num_atom_ord = L

Nomes_lista = []
for item in distancias:
    linha1 = [item[0]] + list(item[2:5])
    Nomes_lista.append(linha1)
    linha1 = [item[1]] + list(item[5:8])
    Nomes_lista.append(linha1)

Nomes = []
for k in range(1, L + 1):
    for item in Nomes_lista:
        if k == int(item[0]):
            Nomes.append(item)
            break
Xo = []

for item in Nomes:
    for atomo in pdb:
        atomo_name = [atomo[2], atomo[3], atomo[5]]
        if all(item[k + 1] == atomo_name[k] for k in range(3)):
            Xo.append(np.array(atomo[6:9], dtype='float'))

# parâmetros referentes ao ponto inicial Xo;
Solution = np.array(Xo)
# Solve :: Solução esperada/correta;
Solve = centralizar(Solution)

# Determinação do ponto inicial:
print(":: Ponto inicial -- responda Sim apenas em caso afirmativo.")
if_relax = input(">> Ponto proveniente de relaxação convexa: ")
gram = []
Xi = []

if if_relax == "Sim":
    print(">> Leitura do arquivo iniciada...")
    try:
        raid1 = "C:\\Users\\viniv\\Desktop\\Testes\\relax\\Teste {}\\relax_{}.txt".format(filename, filename)
        Xi = np.genfromtxt(raid1, dtype='str')
        Xi = centralizar(Xi)
        print(":: Leitura de relax_{} concluida. ".format(filename))
        raid2 = "C:\\Users\\viniv\\Desktop\\Testes\\relax\\Teste {}\\gram_{}.txt".format(filename, filename)
        gram = np.genfromtxt(raid2)
        print(":: Leitura de gram_{} concluida. ".format(filename))
    except:
        print(":: Arquivo relax_{} ou gram_{} não encontrado...\n"
              ":: Gerando ponto inicial através de perturbação da solução esperada.".format(filename, filename))
else:
    Xi = Solve + Noise * np.asarray([np.random.normal(0, 1, len(Solve)) for i in range(3)]).T

Xi = np.array(Xi, dtype='float')
X_input = Xi  # ponto inicial, pode ser solução com ruído ou resultado da relaxação convexa.

# Demais ajustes...

u, v = np.array(distancias[:, 0], dtype='int'), np.array(distancias[:, 1], dtype='int')
for i in range(len(u)):
    u[i] += -1
    v[i] += -1

lb, ub = np.array(distancias[:, 8], dtype='float'), np.array(distancias[:, 9], dtype='float')


# Funções de suporte:

def rmsd(A, B):
    n = len(A)
    A = centralizar(A)
    B = centralizar(B)
    SVD = svd(np.dot(B, A.T))
    Q = np.dot(SVD[0], SVD[2].T)

    correlation = np.linalg.norm(np.dot(Q, A) - B, ord='fro')

    return np.sqrt(1 / n) * correlation


def mde(ponto, vec_u, vec_v):
    dim = len(vec_u)
    S = 0
    for s in range(dim):
        dist_ponto = distance(vec_u[s], vec_v[s], ponto)
        S += max((lb[s] - dist_ponto) / lb[s], 0) + max((dist_ponto - ub[s]) / ub[s], 0)

    return S / m


m = len(lb)
w = np.ones(m)


def dist_matrix_projection(dim, vec_u, vec_v, lower_bound, upper_bound, F):
    y = np.zeros(dim)
    for s in range(dim):
        d = distance(vec_u[s], vec_v[s], F)
        if lower_bound[s] <= d <= upper_bound[s]:
            y[s] = d
        elif d < lower_bound[s]:
            y[s] = lower_bound[s]
        elif d > upper_bound[s]:
            y[s] = upper_bound[s]
    return y


yi = dist_matrix_projection(m, u, v, lb, ub, Xi)
prop_dist = 0
for k in range(m):
    if int(distancias[k][-1]) != 0:
        prop_dist += 1
prop_dist = int(prop_dist)

print("##########################  INFORMAÇÕES  ##########################")
print("Proteína: {}, número de átomos iniciais: {}, após reordenação {}.".format(filename, Num_atom_ini, num_atom_ord))
print('Distâncias avaliadas: {} e Distâncias Conhecidas: {}.'.format(m, prop_dist))
if if_relax == 'Sim':
    print('Ponto inicial proveniente de relaxação convexa.')
    d_gram = dist_matrix_projection(m, u, v, lb, ub, gram)
    fo_gram = stress(gram, d_gram, u, v, w)
    print(':: Valor inicial da função objetivo (problema relaxado): ')
    print('valor obj (Relax): {}'.format(fo_gram))
print(":: Iniciando spg...")
try:
    to = time.time()  # tempo inicial;
    X_spg, BackTracking, iter, fun_o, GtD = spg(stress, grad_stress, Xi, yi, u, v, w, lb, ub, TOL, N, M)
    t = time.time() - to  # duração total to processo.
    X = X_spg
    print("Solução encontrada !")
    print('Iterações: {:<21} valor obj: {}\n'
          'GtD   = {:<24} tempo(s) = {}'.format(iter, fun_o / 2, GtD, t))
except ValueError:
    BackTracking = " -- "
    X_spg = "ERRO"
    X = X_spg
    print("ERRO !!!")

print('RMSDi = {:<24} RMSDf = {}'.format(rmsd(Xi, Solve), rmsd(X, Solve)))
print('MDEi = {:<25} MDEf = {}'.format(mde(Xi, u, v), mde(X, u, v)))
print('BackTracking = {}.'.format(BackTracking))
print("###################################################################")

# Geração do arquivo contendo o ponto inicial utilizado e os arquivos PDB finais:

raid = 'C:\\Users\\viniv\\Desktop\\Testes\\Teste {}\\Ponto_inicial.txt'.format(filename)
with open(raid, 'a') as Ponto:
    for item in X_input:
        item = str(item).replace('[', ' ').replace(']', ' ')
        Ponto.write("%s\n" % item)

# Arquivo de saída PDB:
# ATOM i atom amino A res x1 x2 x3 0.00  0.00 atom[0]
Nomes = np.array(Nomes)


def write_pdb_file(dir, coord, atom, res, amino):
    dir = str(dir)
    with open(dir, 'w') as fid:
        fid.write('HEADER    DE NOVO PROTEIN                         01-JAN-14   XXXX              \n')
        fid.write('TITLE     SOLUTION OBTAINED BY iBP        \n')
        fid.write('REMARK   1  Branch and Prune for the DMDGP        \n')
        fid.write('REMARK   1  \n')
        fid.write('MODEL        1   \n')
        # n: número de átomos;
        n = len(atom)
        # naa: n-ésimo resíduo na sequência de aminoácidos;
        naa = int(res[n - 1])

        if naa > 0:
            if int(res[n - 2]) > naa:
                naa = int(res[n - 2])

            nn = 0

            fid.write('SEQRES   1 A {:.4f}  '.format(naa))
            j = 0
            k = 1

            for i in range(n):
                if int(res[i]) <= nn:
                    continue
                nn = int(res[i])
                j = j + 1
                if j > 13:
                    j = 1
                    fid.write(' \n')
                    k = k + 1
                    fid.write('SEQRES %3d A %4d  ' % (k, naa))

                fid.write('%s' % amino[i])
            fid.write('\n')

        for i in range(n):
            satom = atom[i]
            dres = res[i]
            samino = amino[i]
            fid.write('ATOM  %5d %4s %3s A%4d    %8.3f%8.3f%8.3f  0.00  0.00          %1s \n' % (
                i + 1, satom, samino, int(dres), coord[i, 0], coord[i, 1], coord[i, 2], satom[0]))

        fid.write('TER     %d      %s A  %d \n' % (n, amino[n - 1], int(res[n - 1])))
        fid.write('ENDMDL  \n')
        fid.write('END \n')


raid = 'C:\\Users\\viniv\\Desktop\\Testes\\Teste {}\\Sol_{}.pdb'.format(filename, filename)
write_pdb_file(raid, X_spg, Nomes[:, 1], Nomes[:, 3], Nomes[:, 2])
raid = 'C:\\Users\\viniv\\Desktop\\Testes\\Teste {}\\Orig_{}.pdb'.format(filename, filename)
write_pdb_file(raid, Solve, Nomes[:, 1], Nomes[:, 3], Nomes[:, 2])
Xi = centralizar(Xi)
raid = 'C:\\Users\\viniv\\Desktop\\Testes\\Teste {}\\Ponto_{}.pdb'.format(filename, filename)
write_pdb_file(raid, Xi, Nomes[:, 1], Nomes[:, 3], Nomes[:, 2])
