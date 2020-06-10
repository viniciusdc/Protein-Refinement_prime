import numpy as np
import argparse

parser = argparse.ArgumentParser(description="A simples tool to reformat the Mdjeep distance file format")
parser.add_argument('protein_name', type=str, help="input the Mdjeep protein name")

args = parser.parse_args()

main_dir = f'C:\\Users\\viniv\\Documents\\protein_tests\\Mdjeep'
direct = main_dir + f'\\{args.protein_name}.txt'
distancias = np.genfromtxt(direct, dtype='str')
distancias[:, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]] = distancias[:, [1, 0, 7, 9, 3, 6, 8, 2, 4, 5]]

size = len(distancias[:, 1])
pointer = np.zeros(size, dtype=int)

for i in range(size):
    lb = float(distancias[i, 8])
    ub = float(distancias[i, 9])
    error = ub - lb
    if abs(error) > 1e-2:
        pointer[i] = int(1)
    distancias[i, 8] = str(lb - 1e-6)
    distancias[i, 9] = str(ub + 1e-6)

print(">> Preparação para escrita em arquivo;")
# formatação do arquivo de saída
dist_data = []
k = 0
for item in distancias:
    linha = '  '
    for i in range(len(item)):
        if i in [8, 9]:
            linha += format(item[i]).ljust(23, ' ')
        else:
            linha += format(item[i]).ljust(5, ' ')
    linha = linha + format(str(pointer[k]).ljust(5, ' '))
    dist_data.append(linha)
    k += 1

    # Escrita em arquivo:
output_d = main_dir + f'\\Teste {args.protein_name}\\dist_{args.protein_name}.txt'
with open(output_d, 'w') as f:
    for item in dist_data:
        f.write("%s\n" % item)

