import random as rn
from pprint import pprint
import time

start_time = time.time()


def generate_table(A):
    table = []
    for i in range(8):
        x1 = rn.randint(0, 20)
        x2 = rn.randint(0, 20)
        x3 = rn.randint(0, 20)
        y = A[0] + A[1]*x1 + A[2]*x2 + A[3]*x3

        table.append([x1, x2, x3, y])
    return table

def X0(table):
    x0t = [[], []]
    for i in range(len(table[0]) - 1):
        max = table[0][i]
        min = table[0][i]
        for j in table:
            if j[i] > max:
                max = j[i]
            if j[i] < min:
                min = j[i]
        x0 = (max + min)/2
        Dx = max - x0
        x0t[0].append(x0)
        x0t[1].append(Dx)
    return x0t

def Xn(table, x0t):
    xnt = []
    for i in table:
        line = []
        for j in range(len(table[0]) - 1):
            xn = (i[j] - x0t[0][j])/x0t[1][j]
            line.append(xn)
        xnt.append(line)
    return xnt

def variant_assighment(table, x0t, A):
    ye = A[0] + A[1]*x0t[0][0] + A[2]*x0t[0][1] + A[3]*x0t[0][2]
    print('Ye = {}'.format(ye))
    max = 0
    answer = []
    for i in table:
        res = (i[-1] - ye)**2
        if res > max:
            max = res
            answer = i
    return answer

a0 = 4
a1 = 6
a2 = 9
a3 = 1
A = [a0, a1, a2, a3]

table = generate_table(A)
pprint(table)
print("\n")

X0_table = X0(table)
pprint(X0_table)
print("\n")

Xn_table = Xn(table, X0_table)
pprint(Xn_table)
print("\n")

answer = variant_assighment(table, X0_table, A)
print(answer)


print("час роботи: {}".format(time.time() - start_time))



