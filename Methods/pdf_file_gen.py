"""Função responsável por gerar o arquivo .pdb final da proteina XXX"""


# Arquivo de saída PDB:
# ---------- ATOM   i  atom amino  A  res  x1 x2 x3  --  -- atom[0]
# exemplo :: ATOM   1    N   LYS   A   2   x  y  z   --  --   N

# local_dir :: diretório; caminho onde o arquivo será hospedado;
# coord :: coordenadas dos átomos escolhidos;
# atom :: conjunto de átomos escolhidos;
# res :: resíduos dos átomos escolhidos;
# amino :: aminoacidos dos átomos escolhidos;
def write_pdb_file(local_dir, coord, atom, res, amino):
    local_dir = str(local_dir)
    with open(local_dir, 'w') as fid:
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
            s_atom = atom[i]
            d_res = res[i]
            samino = amino[i]
            fid.write('ATOM  %5d %4s %3s A%4d    %8.3f%8.3f%8.3f  0.00  0.00          %1s \n' % (
                i + 1, s_atom, samino, int(d_res), coord[i, 0], coord[i, 1], coord[i, 2], s_atom[0]))

        fid.write('TER     %d      %s A  %d \n' % (n, amino[n - 1], int(res[n - 1])))
        fid.write('ENDMDL  \n')
        fid.write('END \n')

    return
