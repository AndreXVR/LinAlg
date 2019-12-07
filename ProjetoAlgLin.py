import numpy as np

#INPUT DA MATRIZ
def converteString(string):
  list = string.split(" ")
  list = [int(x) for x in list]
  return list
  
def inputMatriz(m, n):
    A = [] 
    for i in range(m):
        lin = input()
        lin = converteString(lin)
        A.append(lin)
    
    A = np.array(A)
    if np.shape(A) != (m,n):
      return "Dimensoes da matriz invalidas."
    else:
      return A

def matQuad(A):
    m, n = np.shape(A)
    return m == n

#IMPLEMENTACAO DO TRAÇO

def trace(A):
    m, n = np.shape(A)

    if !matQuad(A):
        return "Não é possivel encontrar traço de matrizes não quadradas."
    
    for i in range(m):
        soma = soma + A[i][i]

    return soma

#IMPLEMENTACAO DA DETERMINANTE
def subMatriz(A, lin, col):
    m, n = np.shape(A)
    submatriz = []

    for i in range(m):
        sm_lin = []
        if i != lin:
            for j in range(n):
                if j != col:
                    sm_lin.append(A[i][j])

            submatriz.append(sm_lin)
    submatriz = np.array(submatriz)
    return submatriz 

def determinante(A):  
    det = 0.0
    list_cofat = []
    m, n = np.shape(A)
    if !matQuad(A):
        return "Não é possivel encontrar determinante de matrizes não quadradas."

    i = 0

    if m == 2:
        return (A[0][0]*A[1][1]) - (A[0][1]*A[1][0])

    elif m == 1:
        return A[0][0]

    else:
        for j in range(n):
            list_cofat.append(((-1)**((i+1)+(j+1))) * determinante(subMatriz(A, i, j)))
            det = (list_cofat[j] * A[i][j]) + det

    return det

#IMPLEMENTACAO DA TRANSPOSTA
def transposta(A):
    m, n = np.shape(A)
    tA= np.zeros((n,m),int)
    for i in range(m):
        for j in range(n):
            tA[j][i] = A[i][j]
    return tA

#IMPLEMENTACAO DA INVERSA
def matrizCof(A):
    m, n = np.shape(A)
    cA = np.zeros((m,n),int)
    for i in range(m):
        for j in range(n):
            cA[i][j] = ((-1)**((i+1)+(j+1))) * np.det(subMatriz(A, i, j))

    return cA


def inversa(A):
    if !matQuad(A):
        return "Não é possivel encontrar a inversa de matrizes não quadradas."
    invA = 1/np.det(A) * transposta(matrizCof(A)))
    return invA

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




