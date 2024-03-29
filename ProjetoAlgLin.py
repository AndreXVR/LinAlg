import numpy as np
from SubFunc import *

#INPUT DA MATRIZ

def formalizaNum(num):
    if int(num)!=0:
        if abs(num) < 0.00001:
            return 0
        elif num/int(num) == 1:
            return int(num)
        else:
            return round(num,2)
    else:
        if abs(num) < 0.00001:
            return 0
        else:
            return round(num,2)

def formalizaMat(A):
    m, n = np.shape(A)
    for i in range(m):
        for j in range(n):
            A[i][j] = formalizaNum(A[i][j])
    
    return A

def formalizaVet(V):
    n = len(V)
    for i in range(n):
        V[i] = formalizaNum(V[i])
    
    return V

def converttoNum(string):
    list = []
    for i in range(len(string)):
        list.append(int(string[i]))
    
    return list

  
def inputMatriz(string):
    A = []
    lin = []
    cont = 0
    list = string.split(" ")
    list = converttoNum(list)
    n = len(list)
    m = int(n**(1/2))
    for i in range(m):
        lin = []
        for j in range(m):
            lin.append(list[j+cont])

        cont = cont + m
        A.append(lin)
            
    A = np.array(A)
    return A 

def matQuad(A):
    m, n = np.shape(A)
    return m == n

#IMPLEMENTACAO DO TRAÇO

def trace(A):
    n = np.shape(A)[0]
    soma = 0
    if not matQuad(A):
        return "Não é possivel encontrar traço de matrizes não quadradas."
    
    for i in range(n):
        soma = soma + A[i][i]

    return formalizaNum(soma)

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
    if not matQuad(A):
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

    return formalizaNum(det)

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
            cA[i][j] = ((-1)**((i+1)+(j+1))) * np.linalg.det(subMatriz(A, i, j))

    return cA


def inversa(A):
    if not matQuad(A):
        return "Impossivel com matriz não quadrada."
    if np.linalg.det(A) == 0:
        return "Impossivel com matriz singular"
    invA = 1/np.linalg.det(A) * transposta(matrizCof(A))
    return formalizaMat(invA)

#IMPLEMENTACAO POLINOMIO



def somaDiagP(A, valor):
    m,n = np.shape(A)
    for i in range(m):
        A[i][i] = A[i][i] + valor

    return A

def stringPolinomio(coeficientes,i):
    k = len(coeficientes)
    if i < k-1:
        polinomio = "+" + "(" + str(formalizaNum(coeficientes[i])) + ")" + chr(955) + "^" + str(k-1-i) + " " + stringPolinomio(coeficientes,i+1)
    else:
        polinomio = "+" + "(" + str(formalizaNum(coeficientes[i])) + ")"
    return polinomio

def polinomio(A):
    m, n = np.shape(A)
    if not matQuad(A):
        return "Não é possivel encontrar polinomio caracteristico de matrizes não quadradas."
    coeficientes = np.array([1.])
    auxA = np.array(A)
    for i in range(1, n + 1):
        coef = -auxA.trace() / i
        coeficientes = np.append(coeficientes, coef)
        auxA = somaDiagP(auxA, coef)
        auxA = np.dot(A, auxA)
        
    return stringPolinomio(coeficientes,0)

#IMPLEMENTACAO AUTOVALORES

def fatorizacaoQR(A):
    produtoQR = A
    m, n = np.shape(produtoQR)
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
    m, n = np.shape(A)
    strRetorno = ""
    if not matQuad(A):
        return "Não é possivel encontrar autovalores de matrizes não quadradas."
    listAutovalor = []
    for i in range(m):
        autovalor = A[i][i]
        formalizaNum(autovalor)
        listAutovalor.append(formalizaNum(av(A)[i]))
        strRetorno = strRetorno +chr(955)+" = " + str(listAutovalor[i]) + "\n"

    return strRetorno

#IMPLEMENTACAO AUTOVETORES

def calcB(A):
    m, n = np.shape(A)
    auxA = A
    B = []
    for i in range(1, n):
        trace = -auxA.trace() / i 
        auxA = somaDiagP(auxA, trace)
        B.append(auxA)
        auxA = np.dot(A, auxA) 
    return B

def kColuna(A, k):
    m, n = np.shape(A)
    col = []
    for i in range(n):
        num = []
        num.append(A[i][k])
        col.append(num)
    return col

def verifVet(vet):
    n = len(vet)
    for i in range(n):
        if abs(vet[i]) > 0.000001:
            return 1
    return 0
    
def diagonal(A):
  m, n = np.shape(A)
  for i in range(m):
    for j in range(n):
      if i != j:
        if A[i][j] != 0:
          return False
  
  return True
  

def autovetores(A):
    strRetorno = ""
    if not matQuad(A):
        return "Não é possivel encontrar autovetores de matrizes não quadradas."
    avals = np.linalg.eigvals(A)
    B = calcB(A)
    m, n = np.shape(A)
    id = np.identity(m, int)
    numAvals = len(avals)
    numBs = len(B)
    cont = 0
    k = 0
    listAutovetores = []
    for i in range(m):
        listAutovetores.append((formalizaNum(ave(A)[0][i]),formalizaVet(ave(A)[1][i])))
        strRetorno = strRetorno +chr(955)+'='+str(listAutovetores[i][0]) +" "+str(listAutovetores[i][1])+'\n'
    while cont != numAvals:
        u = np.array(kColuna(id, k))
        for j in range(numBs):
            u = ((avals[cont]*u) + np.array(kColuna(B[j], k)))
        if verifVet(u) == 0:
            k+=1
            u = np.array(kColuna(id, k))
        else:
            cont+=1
            k = 0
    
    return strRetorno

def matrizD(A):
    if not matQuad(A):
        return "Não é possivel encontrar matriz diagonal de matrizes não quadradas."

    n = np.shape(A)[0]
    if len(np.linalg.eigvals(A)) != n:
        return "A matriz é defectiva"

    D = np.zeros((n,n),int)
    for i in range(n):
        D[i][i] = formalizaNum(np.linalg.eigvals(A)[i])
        
    return D
            


