import numpy as np

def le_mat(m, n):
    mat = []
    for i in range(m):
        mat.append([])
        for j in range(n):
            print("Valor na pos ", i + 1, ",", j + 1)
            num = float(input())
            mat[i].append(num)
    return mat


def mat_escalar(escalar, mat):
    mat_prod = []
    n_lin = len(mat)
    n_col = len(mat[0])
    for i in range(n_lin):
        mat_prod.append([])
        for j in range(n_col):
            mat_prod[i].append(mat[i][j]*escalar)
    return mat_prod


def mat_menor(mat, lin, col):  # Recebe uma matriz A e uma posiçao da matriz, retorna a determinante da matriz formada por A menos a linha i e a coluna j.
    menor = []  # cria uma matriz
    n_lin = len(mat)  # diz a quantidade de linhas
    n_col = n_lin  # Quantidade de colunas

    for i in range(n_lin):
        list_lin = []

        for j in range(n_col):
            if j != col:
                list_lin.append(mat[i][j])

        if i != lin:
            menor.append(list_lin)

    return menor


def mat_det(mat):  # Recebe uma matriz A e uma posiçao da matriz, retorna a determinante da matriz formada por A menos a linha i e a coluna j.
    det_val= 0.0
    list_cofat = []
    n_lin = len(mat)
    n_col = n_lin
    i = 0

    if n_lin == 2:
        return (mat[0][0]*mat[1][1]) - (mat[0][1]*mat[1][0])

    elif n_lin == 1:
        return mat[0][0]

    else:
        for j in range(n_col):
            list_cofat.append(((-1)**((i+1)+(j+1))) * mat_det(mat_menor(mat, i, j)))
            det_val = (list_cofat[j] * mat[i][j]) + det_val

    return det_val


def mat_transposta(mat):
    transposta=[]
    n_lin= len(mat)
    n_col= len(mat[0])
    for i in range(n_col):
        transposta.append([])
        for j in range(n_lin):
            transposta[i].append(mat[j][i])
    return transposta


def mat_cofatores(mat):
    comat = []
    n_lin = len(mat)
    n_col = n_lin

    for i in range(n_lin):
        comat.append([])

        for j in range(n_col):
            comat[i].append(((-1)**((i+1)+(j+1))) * mat_det(mat_menor(mat, i, j)))

    return comat


def mat_inversa(mat):
    inversa = mat_escalar(1/mat_det(mat), mat_transposta(mat_cofatores(mat)))
    return inversa

#IMPLEMENTACAO POLINOMIO

def polinomio(A): '''calcula usando faddeev leverrier'''
    n = len(A[0])
    A = np.array(A)
    coeficientes = np.array([1.])
    A2 = np.array(A)
    for i in range(1, n + 1):
        coef = -A2.trace() / i
        coeficientes = np.append(coeficientes, coef)
        A2 = subtraiDiagP(A2, coef)
        A2 = np.dot(A, A2)
    return coeficientes

def subtraiDiagP(mat, num):
    n_lin = len(mat)
    for i in range(n_lin):
        mat[i][i] = mat[i][i] + num

    return mat '''retorna uma lista de autovalores, onde a posicao 0 é x^n, sendo n a ordem da matriz, a posicão 1, x^n-1 e assim por diante'''

#IMPLEMENTACAO AUTOVALORES

import numpy as np


def fatorizacaoQR(A):
    produtoQR = np.array(A)
    m, n = produtoQR.shape
    it = 1000
    k = 0
    while k < it:
        k = k + 1
        Q = np.zeros((m, n))
        R = np.zeros((n, n))
        produtoQR = np.transpose(produtoQR)
        R[0][0] = np.linalg.norm(produtoQR[0], 2)
        Q[:, 0] = produtoQR[0] / R[0][0]
        for i in range(1, n):
            Q[:, i] = produtoQR[i]
            for j in range(0, i):
                R[j][i] = np.dot(Q[:, j], Q[:, i])
                Q[:, i] = Q[:, i] - (R[j][i] * Q[:, j])
            R[i][i] = np.linalg.norm(Q[:, i], 2)
            Q[:, i] = Q[:, i] / R[i][i]

        produtoQR = np.dot(R, Q)

    return produtoQR


def autovalores(A):
    n = A.shape[0]
    listAutovalor = []
    for i in range(n):
        autovalor = A[i][i]
        if (abs(autovalor) < 0.00000001):  # Faz com que autovalores aproximados de zero recebam 0
            autovalor = 0
        listAutovalor.append(autovalor)
    return listAutovalor

#IMPLEMENTACAO AUTOVETORES

import numpy as np


def subtraiDiagP(mat, num):
    n_lin = len(mat)
    for i in range(n_lin):
        mat[i][i] = mat[i][i] + num
    return mat


def matrizBs(mat):
    n = len(mat[0])
    mat = np.array(mat)
    mat2 = np.array(mat)
    matB = []
    for i in range(1, n):
        trac = -mat2.trace() / i
        mat2 = subtraiDiagP(mat2, trac)
        matB.append(mat2)
        mat2 = np.dot(mat, mat2)
    return matB


def colunaDeMat(mat, posCol):
    n = len(mat)
    col = []
    for i in range(n):
        a = []
        a.append(mat[i][posCol])
        col.append(a)
    return col


def verifVet(vet):
    n = len(vet)
    for i in range(n):
        if vet[i] != 0:
            return 1
    return 0


def autovetor(mat):
    autv = np.linalg.eigvals(mat)
    matB = matrizBs(mat)
    tam_mat = len(mat)
    id = np.identity(tam_mat, int)
    qtd = len(autv)
    qtd_matB = len(matB)
    cont = 0
    col = 0
    result = []
    while cont != qtd:
        u = np.array(colunaDeMat(id, col))
        for j in range(qtd_matB):
            u = ((autv[cont] * u) + np.array(colunaDeMat(matB[j], col)))
        if verifVet(u) == 0:
            col += 1
            u = np.array(colunaDeMat(id, col))
        else:
            result.append(u)
            cont += 1
            col = 0
    return result


mat = [[2, -2, 3], [0, 3, -2], [0, -1, 2]]
m = autovetor(mat)
print(m)'''




