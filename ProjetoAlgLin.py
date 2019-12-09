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

def formalizaNum(num):
    if abs(num) < 0.0000001:
        return 0
    elif num/int(num) == 1:
        return int(num)
    else:
        return round(num,2)

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
    if !matQuad(A):
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
    m, n = A.shape
    if !matQuad(A):
        return "Não é possivel encontrar autovalores de matrizes não quadradas."
    listAutovalor = []
    for i in range(m):
        autovalor = A[i][i]
        formalizaNum(autovalor)
        listAutovalor.append(autovalor)
    return listAutovalor

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
    if !matQuad(A):
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
    if diagonal(A):
      for i in range(m):
        listAutovetores.append(id[i][:])
    else:
      while cont != numAvals:
          u = np.array(kColuna(id, k))
          for j in range(numBs):
              u = ((avals[cont]*u) + np.array(kColuna(B[j], k)))
          if verifVet(u) == 0:
              k+=1
              u = np.array(kColuna(id, k))
          else:
              listAutovetores.append(u)
              cont+=1
              k = 0
    
    return listAutovetores
  
  #Codigo Interface
from numpy import *
from Funcoes import *
from functools import partial
from tkinter import *

window = Tk()
window["bg"] = "black"
window.title("Operações com Matrizes")
window.geometry("1300x500+0+0")
lb2 = Label(window, text="Insira a matriz em python",
font="arial 20 bold", width=45, bg="orange", fg="white")
lb2.grid(row = 0, column = 0)
lbaux = Label(window, font="arial 1 ", bg = "black")
lbaux.grid(row = 1, column = 0)
edmatriz = Entry(window, width = 30)
edmatriz.insert(INSERT, "[[0,0,0],"
"[0,0,0],"
"[0,0,0]]")
edmatriz.grid(row=2, column = 0)

lbaux2 = Label(window, font="arial 1 bold", bg = "black")
lbaux2.grid(row = 3, column = 0)
bt = Button(window,width=19,text="Determinante", relief=GROOVE,font=" arial 12", bg =
"orange", activebackground="green", fg="white", command = determinante)
bt.grid(row=4,column = 0)
lbaux3 = Label(window, font="arial 1 bold", bg = "black")
lbaux3.grid(row = 5, column = 0)
bt2 = Button(window,width=19,text="Traço", relief=GROOVE, font=" arial 12", bg = 
"orange", activebackground="green", fg="white")
#bt2.grid(row=12,column=0)
bt2.place(x=450,y=85)
lbaux4 = Label(window, font="arial 1 bold", bg = "black")
lbaux4.grid(row = 7, column = 0)
bt3 = Button(window,width=19,text="Transposta", relief=GROOVE, font=" arial 12", bg =
"orange", activebackground="green", fg="white")
bt3.grid(row=6,column = 0)
lbaux5 = Label(window, font="arial 1 bold", bg = "black")
lbaux5.grid(row = 9, column = 0)
bt4 = Button(window,width=19,text="Inversa", relief=GROOVE, font=" arial 12", bg =
"orange", activebackground="green", fg="white")
#bt4.grid(row=14,column=0)
bt4.place(x=450,y=121)
lbaux6 = Label(window, font="arial 1 bold", bg = "black")
lbaux6.grid(row = 11, column = 0)
bt5 = Button(window,width=19,text="Polinômio Característico", relief=GROOVE, font=" arial 12",
bg = "orange", activebackground="green", fg="white")
bt5.grid(row=8,column = 0)
lbaux7 = Label(window, font="arial 1 bold", bg = "black")
lbaux7.grid(row = 13, column = 0)
bt6 = Button(window,width=19,text="Autovalores", relief=GROOVE, font=" arial 12", bg =
"orange", activebackground="green", fg="white")
#bt6.grid(row=16,column=0)
bt6.place(x=450,y=157)
lbaux8 = Label(window, font="arial 1 bold", bg = "black")
lbaux8.grid(row = 15, column = 0)
bt7 = Button(window,width=19,text="Autovetores", relief=GROOVE, font=" arial 12", bg =
"orange", activebackground="green", fg="white")
bt7.grid(row=10,column = 0)
lbaux9 = Label(window, font="arial 1 bold", bg = "black")
lbaux9.grid(row = 17, column = 0)
bt8 = Button(window,width=19,text="Matriz Diagonal", relief=GROOVE, font=" arial 12", bg =
"orange", activebackground="green", fg="white")
#bt8.grid(row=18,column=0)
bt8.place(x=450,y=193)
lbaux10 = Label(window, font="arial 1 bold", bg = "black")
lbaux10.grid(row = 19, column = 0)

lb2_0 = Label(window, text="Resultado", font="arial 20 bold", width=45, bg="orange",
fg="white")
lb2_0.grid(row = 0, column = 2)
lbaux13 = Label(window, font="arial 1 ", bg = "black")
lbaux13.grid(row = 1, column = 2)
ed2 = Entry(window, width = 30)
ed2.insert(INSERT, "Resultado da Operação: ")
ed2.grid(row=3, column = 2)
#lbaux14 = Label(window, font="arial 1 ", bg = "black")
#lbaux14.grid(row = 3, column = 2)
#lbaux15 = Label(window, font="arial 12 bold", width=24, bg="orange", fg="white",
#text="Resultado da Operação:")
#lbaux15.grid(row = 4, column = 2)
#lbaux16 = Label(window, font="arial 1 ", bg = "black")
#lbaux16.grid(row = 5, column = 2)
#lbaux17 = Label(window, font="arial 12 bold ", bg = "white",fg="white", width=40)
#lbaux17.grid(row = 6, column = 2)

window.maxsize(width=1520, height=800)
window.mainloop()
