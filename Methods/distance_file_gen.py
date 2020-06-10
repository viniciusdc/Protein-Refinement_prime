import numpy as np

"O objetivo deste modulo é armazenar a função que cria e gerencia os arquivos de distâncias provenientes do arquivo pdb"


# pdb :: arquivo contendo nome, cadeia e posição dos átomos da proteina XXX
# exemplo :: ['ATOM', '1', 'N', 'LYS', 'A', 'x', 'y', 'z', '-', '-', 'N']
# protein_id :: nome da proteina;
# raid :: diretório local onde o arquivo de saída será hospedado (str).
# exemplo :: 'C:\\Users\\**\\Desktop\\Testes'
def gen_distance_file(pdb, protein_id, raid, overwrite=False):
    # Informações de controle iniciais ---------------------------

    # raid_d :: diretório de saída do arquivo de distâncias gerado
    raid_d = raid + '\\Teste {}\\dist_{}.txt'.format(protein_id, protein_id)

    try:
        with open(raid_d) as f:
            # Se o arquivo de distâncias já existe não há necessidade em criar um novo;
            if overwrite:
                # Usuario requisitou a geração de um novo arquivo de distâncias;
                pass
            else:
                return str(raid_d)

    # caso o arquivo não seja encontrado e/ou não exista no diretório especificado, prosseguimos criando um novo;
    except FileNotFoundError:
        print("Arquivo não encontrado. Gerando novo arquivo !")

    # delta :: ruído adicionado aos intervalos de distâncias exatas, evitando problemas de DL.
    # distance_accept_control :: distância intermolecular limite (Ang) aceita;
    delta = np.random.rand(1)[0] / 10000
    distance_accept_control = 5.0

    print(":: Grau de aceitação da Distância intermolecular: {} Ang.".format(distance_accept_control))

    # Inicio do código -------------------------------------------

    # atomos_aceitos :: lista de átomos aceitos para comparação;
    # cadeia_principal :: definindo o que são os átomos da cadeia principal (Backbone);
    atomos_aceitos = ['H', 'C', 'CA', 'N', 'HA', 'CB', 'O', 'H1']
    cadeia_principal = ['C', 'CA', 'N']

    print('>> Iniciando geração do arquivo de distâncias.')

    data_file = []
    for atom1 in pdb[1:]:
        for atom2 in pdb:
            # Para cada par de átomos (atom1, atom2) verificamos se são átomos aceitáveis;
            if atom1[2] and atom2[2] in atomos_aceitos:
                # Iremos apenas adicionar as distâncias entre átomos diferentes entre sí, e além disso
                # para evitar de adicionar repetidamente a mesma distância, faremos pares do menor para maior:
                if atom1[1] != atom2[1] and int(atom2[1]) < int(atom1[1]):
                    # Cálculo da distância entre os dois átomos:

                    # zyz1 :: posições do átomo correspondente a linha de atom1;
                    # xyz2 :: posições do átomo correspondente a linha de atom2;
                    xyz1 = np.array(atom1[6:9], dtype=float)
                    xyz2 = np.array(atom2[6:9], dtype=float)

                    # dist :: distância (Euclidiana) inicial entre os átomos;
                    dist = np.linalg.norm(xyz1 - xyz2)

                    # Verificaremos se são átomos da cadeia principal.
                    if atom1[2] and atom2[2] in cadeia_principal:
                        # Se estiverem na mesma cadeia -- [1]:
                        if atom1[5] == atom2[5]:
                            # adicionamos um ruído delta na distância para evitar LD;
                            dist_lb = dist - delta
                            dist_ub = dist + delta
                            # por fim, adicionamos esse par e essas distâncias ao arquivo de distâncias:
                            # exemplo : ['1', '6', 'N', 'LYS', '2', 'N', 'LEU', '3', 'dist_l', 'dist_u', '1'];
                            data_file.append([atom2[1], atom1[1], atom2[2], atom2[3],
                                              atom2[5], atom1[2], atom1[3], atom1[5],
                                              dist_lb, dist_ub, '1'])
                        # Caso contrário, adicionaremos esse par se "dist" entre eles for menor que o aceitável -- [0]:
                        elif dist <= distance_accept_control:
                            # Adicionando 'ruído' as limitações inferiores e superiores:
                            dist_lb = dist - 1 / 3
                            dist_ub = dist + 1 / 3
                            # por fim, adicionamos esse par e essas distâncias ao arquivo de distâncias:
                            # exemplo : ['1', '6', 'N', 'LYS', '2', 'N', 'LEU', '3', 'dist_l', 'dist_u', '0'];
                            data_file.append([atom2[1], atom1[1], atom2[2], atom2[3],
                                              atom2[5], atom1[2], atom1[3], atom1[5],
                                              dist_lb, dist_ub, '0'])

    print(">> Reorganizando numeração dos átomos escolhidos;")

    # nesta etapa iremos renumerar os indices dos átomos de acordo com a nova quantidade, pois
    # nem todos os átomos foram aceitos na etapa anterior, logo há necessidade de reindexarmos.

    # index_set :: conjunto dos indices que serão re-indexados
    index_set = np.array(data_file)

    # como o arquivo data_file possui duas colunas de indices, concatenaremos ambas para facilitar o processo;
    index_set = np.array(list(set(list(index_set[:, 0]) + list(index_set[:, 1]))), dtype='int')
    index_set.sort()

    # Agora que temos a lista de todos os indices dos átomos aceitos,sabemos quantos e quais deles foram aceitos;
    # como eventualmente essa lista não conterá todos os átomos da proteína escolhida, haverá "saltos"
    # nos indices listados, portanto é necessária a reordenação desses, de modo a resolver esse problema.

    # O que faremos a seguir é literalmente criar uma nova sequência de indices em data_file, totalmente ordenados;
    # usaremos um dicionário para mapear a nova ordenação, exemplo : {(new_index, atom_index)}
    dict_index_order = enumerate(index_set)

    print(">> Preparação para escrita em arquivo;")
    # formatação do arquivo de saída
    dist_data = []
    for item in data_file:
        linha = '  '
        for i in range(len(item)):
            if i in [8, 9]:
                linha += format(item[i]).ljust(23, ' ')
            else:
                linha += format(item[i]).ljust(5, ' ')
        dist_data.append(linha)

    # Escrita em arquivo:
    with open(raid_d, 'w') as distancias:
        for item in dist_data:
            distancias.write("%s\n" % item)

    # Fim do processo -------------------------------------------------------------------------
    # uma vantagem de criar/ recriar o arquivo de distâncias é o dicionario com o mapa da nova ordenação;
    if overwrite == 1:
        return str(raid_d), dict_index_order
    else:
        return str(raid_d)
