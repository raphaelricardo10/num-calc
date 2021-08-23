import numpy as np
import csv

def seidel_method(A, x, b, err=0.001, prec=3, iteration=0):
    #Calcula o tamanho do sistema
    n = len(A.shape) + 1

    #Monta a matriz trangular superior
    U = np.zeros(shape=(n,n))
    for i in range(0, n):
        for j in range(0, n):
            if j > i:
                U[i][j] = A[i][j]

    #Monta a matriz triangular inferior
    L = A - U

    #Isola as variaveis
    Linv = np.linalg.inv(L)
    Ux = U.dot(x)

    #Resultado = L^-1*(b - Ux)
    xf = Linv.dot(b - Ux)

    #Arredonda o resultado
    xf = xf.round(prec)

    #CSV para análise dos resultados
    with open('seidel.csv', 'a', newline='\n') as csvfile:
        writter = csv.writer(csvfile)
        cols = [iteration, xf[0][0], xf[1][0], xf[2][0]]
        writter.writerow(cols)


    print("Valor de Xf")
    print(xf)

    #Análise de convergência
    conv = True
    for i in range(0, n):
        if round(xf[i][0], prec) != round(x[i][0], prec):
            conv = False

    if conv:
        return xf
    else:
        #Caso não convirja, é feita a análise de erro
        Ax = A.dot(x)
        ok = True
        for i in range(0, n):
            if Ax[i] > b[i]+b[i]*err or Ax[i] < b[i]-b[i]*err:
                ok = False
        if not ok:
        #Realiza outra iteração
            return seidel_method(A, xf, b, iteration=iteration+1)
        else:
        #Erro aceitável
            return xf

#Matriz de coeficientes
A = np.array([[3, -0.1, -0.2],
               [0.1, 7, -0.3],
               [0.3, -0.2, 10]])

#Matriz de variaveis
x = np.array([[0],  #x1
              [0],  #x2
              [0]]) #x3

#Matriz de constantes
b = np.array([[7.85],
              [-19.3],
              [71.4]])

resposta = seidel_method(A, x, b)
