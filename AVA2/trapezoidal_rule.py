import numpy as np
import sympy as sp
import math

def simple_method(f, a, b, var):
    #Cria uma função para que seja avaliada em um ponto
    primitive = sp.lambdify(var, f)


    #Calcula a largura do trapézio
    dx = b-a
    #Reliza a integração
    area = dx*(primitive(a) + primitive(b))/2

    #Encontra f''(x)
    f_2prime = sp.diff(sp.diff(f))
    #Torna f''(x) aplicável num ponto
    f_2prime_P = sp.lambdify(var, f_2prime)

    #Encontra o máximo de f''(csi)
    maxx = f_2prime_P(a)
    for i in range(a,b+1):
        if f_2prime_P(i) > maxx:
            maxx = f_2prime_P(i)

    #Calcula o erro
    error = dx**3*maxx/12
    return [area, error]

def rep_method(f, a, b, var, n):
    #Cria uma função para que seja avaliada em um ponto
    primitive = sp.lambdify(var, f)

    #Calcula a largura de cada trapézio
    dx = (b-a)/n

    #Realiza a integração
    acc = 0
    for k in range(1, n):
        acc += primitive(a + k*dx)
    area = dx*(acc + (primitive(b) + primitive(a))/2)

    #Encontra f''(x)
    f_2prime = sp.diff(sp.diff(f))
    #Torna f''(x) aplicável num ponto
    f_2prime_P = sp.lambdify(var, f_2prime)

    #Encontra o máximo de f''(csi)
    maxx = f_2prime_P(a)
    for i in range(a,b+1):
        if f_2prime_P(i) > maxx:
            maxx = f_2prime_P(i)

    #Calcula o erro
    error = dx**3*maxx/12
    return [area, error]


if __name__ == "__main__":
    #Inicializa a função e a variável
    x = sp.Symbol('x')
    f = sp.sqrt(x)

    #Inicializa os limites de integração
    a = 1
    b = 4

    area1, error1 = simple_method(f, a, b, x)
    area2, error2 = rep_method(f, a, b, x, 5)

    print("Área no método simples: ", area1)
    print("O erro no método simples é de: ", error1)
    print("Área no metodo repetido: ", area2)
    print("O erro no método repetido é de: ", error2)
    input()