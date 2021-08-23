import numpy as np
import csv

def jacobi_method(A, x, b, err=0.001, prec=3, iteration=0, max=1000, dominant=False):
    #Caso tenha alcançado o limite de iterações
    if iteration+1 == max:
        return x
 
    #Calcula o tamanho do sistema
    n = len(A.shape) + 1

    if iteration == 0:
        #Força os tamanhos corretos das matrizes
        A = A.reshape(n,n)
        x = x.reshape(n,)
        b = b.reshape(n,)

        #Checa se a matrix é diagonalmente dominante
        dominant = True
        for i in range(0, n):
            acc = 0
            for j in range(0,n):
                if i != j:
                    acc += abs(A[i][j])
            if A[i][i] < acc:
                dominant = False
                break

    #Monta a matriz diagonal
    lst = []
    for i in range(0, n):
        lst.append(A[i][i])
    D = np.diag(lst)

    #Isola as variaveis
    R = A - D
    Dinv = np.linalg.inv(D)
    Rx = R.dot(x)

    #Resultado = D^-1*(b - Rx)
    xf = Dinv.dot(b-Rx)

    #Arredonda o resultado
    xf = xf.round(prec)
    x = x.round(prec)

    #CSV para análise dos resultados
    with open('jacobi.csv', 'a', newline='\n') as csvfile:
        writter = csv.writer(csvfile)
        cols = [iteration, xf[0], xf[1], xf[2]]
        writter.writerow(cols)

    #print("Valor de Xf")
    #print(xf)

    #Se a matriz for dominante, o critério de parada será o valor exato
    if dominant:
        Ax = A.dot(xf).round(prec)
        for i in range(0, n):
            if Ax[i] != b[i]:
                #Caso o resultado não seja encontrado, realiza outra iteração
                return jacobi_method(A, xf, b, iteration=iteration+1, dominant=True)
            #Retorna o resultado exato
            return xf
    else:
        #Caso não seja, o critério de parada será de convergência ou erro aceitável
        #Análise de convergência
        for i in range(0, n):
            if abs(xf[i]) > abs(x[i]+x[i]*err) or abs(xf[i]) < abs(x[i]-x[i]*err):
                conv = False
                break
        if conv:
            return xf
        else:
            #Caso não convirja, é feita a análise de erro
            Ax = A.dot(xf).round(prec)
            ok = True
            for i in range(0, n):
                if abs(Ax[i]) > abs(b[i]+b[i]*err) or abs(Ax[i]) < abs(b[i]-b[i]*err):
                    ok = False
                    break
            if not ok:
                #Realiza outra iteração
                return jacobi_method(A, xf, b, iteration=iteration+1)
            else:
            #Erro aceitável
                return xf

if __name__ == "__main__" :
    #Matriz de coeficientes
    A = np.array([[3, -0.1, -0.2],
                  [0.1, 7, -0.3],
                  [0.3, -0.2, 10]], np.float32)

    #Matriz de variaveis
    x = np.array([[0],              #x1
                  [0],              #x2
                  [0]], np.float32) #x3

    #Matriz de constantes
    b = np.array([[7.85],
                  [-19.3],
                  [71.4]], np.float32)

    print("Bem vindo!\nA matrizes padrão já foram definidas.\nSão elas:")
    print("A:\n", A)
    print("x:\n", x)
    print("b:\n", b)

    resp = input("Gostaria de alterá-las?\nInsira 0 para sim\nAperte qualquer outra tecla para nao\n")

    if resp == '0':
        n = int(input("Insira a quantidade de variáveis do sistema: "))

        #Letra a
        coeff = 97
        lstA = []
        lstB = []
        lstX = []

        print("A seguir serão solicitados os valores os coeficientes de cada equação do sistema.")
        for i in range(0,n):
            tempLst = []
            for j in range(0, n):
                tempLst.append(float(input(str.format("Insira o valor de {0}: ", chr(coeff+j)))))
            lstB.append(float(input("Insira o valor do termo idependente: ")))
            lstA.append(tempLst)

        for i in range(0, n):
            lstX.append(float(input(str.format("Insira o valor inicial de x{0}: ", i+1))))

        A = np.asarray(lstA, dtype=np.float32)
        x = np.array(lstX, dtype=np.float32)
        b = np.array(lstB, dtype=np.float32)


    resposta = jacobi_method(A, x, b)
    print(resposta)