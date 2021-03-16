import random as rn
from math import sqrt
from numpy import linalg as lg
import pprint

def experiment(x1, x2, y):
    y_gen = [rn.randint(y[0], y[1]) for i in range(5)]
    y_mid = sum(y_gen)/len(y_gen)
    sigma = 0
    for i in range(len(y_gen)):
        sigma += ((y_gen[i] - y_mid)**2)/5
    return [x1, x2, y_gen, y_mid, sigma]

def kriteria(x1, x2, y):
    flag = 1
    m = 5
    Sigma0 = sqrt((2*(2*m - 2))/(m*(m-4)))
    while flag:

        exp1 = experiment(x1[0], x2[0], y)
        exp2 = experiment(x1[1], x2[0], y)
        exp3 = experiment(x1[0], x2[1], y)

        F1 = exp1[-1]/exp2[-1]
        F2 = exp3[-1]/exp1[-1]
        F3 = exp3[-1]/exp2[-1]

        F = [F1, F2, F3]

        o1 = ((m-2)/m)*F1
        o2 = ((m-2)/m)*F2
        o3 = ((m-2)/m)*F3

        O = [o1, o2, o3]

        R1 = abs(o1 - 1)/Sigma0
        R2 = abs(o2 - 1)/Sigma0
        R3 = abs(o3 - 1)/Sigma0

        R = [R1, R2, R3]

        if R1 < 2 and R2 < 2 and R3 < 2:
            flag = 0

    return [[exp1, exp2, exp3], F, O, R]

def koeficients(table):
    mx1 = (-1+1+(-1))/3
    mx2 = (-1+(-1)+1)/3
    my = (table[0][3] + table[1][3] + table[2][3])/3
    a1 = (1+1+1)/3
    a2 = (1-1-1)/3
    a3 = (1+1+1)/3
    a11 = ((-1)*table[0][3] + 1*table[1][3] + (-1)*table[2][3])/3
    a22 = ((-1)*table[0][3] + (-1)*table[1][3] + 1*table[2][3])/3

    b0 = (lg.det([[my, mx1, mx2],
                [a11, a1, a2],
                [a22, a2, a3]]))/(lg.det([[1, mx1, mx2],
                                                [mx1, a1, a2],
                                                [mx2, a2, a3]]))
    b1 = (lg.det([[1, my, mx2],
                  [mx1, a11, a2],
                  [mx2, a22, a3]]))/(lg.det([[1, mx1, mx2],
                                             [mx1, a1, a2],
                                             [mx2, a2, a3]]))
    b2 = (lg.det([[1, mx1, my],
                  [mx1, a1, a11],
                  [mx2, a2, a22]]))/(lg.det([[1, mx1, mx2],
                                             [mx1, a1, a2],
                                             [mx2, a2, a3]]))
    return [b0, b1, b2]

def check(koeficients, x1 = [-1, 1], x2 = [-1, 1]):
    check1 = koeficients[0] + x1[0]*koeficients[1] + x2[0]*koeficients[2]
    check2 = koeficients[0] + x1[1]*koeficients[1] + x2[0]*koeficients[2]
    check3 = koeficients[0] + x1[0]*koeficients[1] + x2[1]*koeficients[2]

    return [ check1, check2, check3]

def naturalisation(x1, x2, koeficients):
    dx1 = abs(x1[1] - x1[0])/2
    dx2 = abs(x2[1] - x2[0])/2
    x10 = (x1[1] + x1[0])/2
    x20 = (x2[1] + x2[0])/2

    a0 = koeficients[0] - koeficients[1]*x10/dx1 - koeficients[2]*x20/dx2
    a1 = koeficients[1]/dx1
    a2 = koeficients[2]/dx2

    return [a0, a1, a2]


x1 = [-30, 20]
x2 = [-70, -10]
print("x1min, x1max = {0}, {1}\n".format(x1[0], x1[1]))
print("x2min, x2max = {0}, {1}\n".format(x2[0], x2[1]))

ymax = (30 - 120)*10
ymin = (20 - 120)*10
y = [ymin, ymax]
print("ymin, ymax = {0}, {1}\n".format(y[0], y[1]))

research = kriteria(x1, x2, y)
print("Fuv:")
print(research[1])
print(" ")

print("Ouv:")
print(research[2])
print(" ")

print("Ruv:")
print(research[3])
print(" ")





table = research[0]
koefs = koeficients(table)
natural = naturalisation(x1, x2, koefs)
check1 = check(koefs)
check2 = check(natural, x1, x2)

print("Таблиця експерименту(x12, yn, ymid, квадратичне відхилення)")
pprint.pprint(table)
print(" ")

print("b0, b1, b2")
print(koefs)
print(" ")

print("перевірка b")
print(check1)
print(" ")

print("перевірка а")
print(check2)
print(" ")
