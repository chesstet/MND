import random as rn
from math import sqrt
import pprint

def exp_row(x1_num, x2_num, x3_num, y, m = 3):
    y_gen = [y[0] + (y[1]-y[0])*rn.random() for i in range(m)]
    y_mid = sum(y_gen)/m
    sigma = 0
    for i in range(len(y_gen)):
        sigma += ((y_gen[i] - y_mid)**2)/m
    return [x1_num, x2_num, x3_num, x1_num*x2_num, x1_num*x3_num, x2_num*x3_num, x1_num*x2_num*x3_num, y_gen, y_mid, sigma]

def experiment(x1, x2, x3, y, m = 3, N = 8):
    flag = 1
    counter = 1
    while flag:
        flag = 0
        exp1 = exp_row(x1[0], x2[0], x3[0], y)
        exp2 = exp_row(x1[0], x2[0], x3[1], y)
        exp3 = exp_row(x1[0], x2[1], x3[0], y)
        exp4 = exp_row(x1[0], x2[1], x3[1], y)
        exp5 = exp_row(x1[1], x2[0], x3[0], y)
        exp6 = exp_row(x1[1], x2[0], x3[1], y)
        exp7 = exp_row(x1[1], x2[1], x3[0], y)
        exp8 = exp_row(x1[1], x2[1], x3[1], y)

        table = [exp1, exp2, exp3, exp4, exp5, exp6, exp7, exp8]
        norm_table = normalize_table(table, x1, x2, x3)


        cochrane = cochrane_kriteria(table, 8)
        cochrane_check = cochrane[0]
        Gp = cochrane[1]

        if not cochrane_check:
            flag = 1
            continue
        else:
            b0 = 0
            b1 = 0
            b2 = 0
            b3 = 0
            b12 = 0
            b13 = 0
            b23 = 0
            b123 = 0

            for i in range(N):
                b0 += 1/N*(norm_table[i][-2])
                b1 += 1/N*(norm_table[i][-2]*norm_table[i][0])
                b2 += 1/N*(norm_table[i][-2]*norm_table[i][1])
                b3 += 1/N*(norm_table[i][-2]*norm_table[i][2])
                b12 += 1/N*(norm_table[i][-2]*norm_table[i][3])
                b13 += 1/N*(norm_table[i][-2]*norm_table[i][4])
                b23 += 1/N*(norm_table[i][-2]*norm_table[i][5])
                b123 += 1/N*(norm_table[i][-2]*norm_table[i][6])

            b = [b0, b1, b2, b3, b12, b13, b23, b123]

            check_b = check(b, table)

            student = student_kriteria(table, 8)
            indexes = student[0]
            SB = student[1]
            Sbeta = student[2]
            student_b = list(map(lambda x: x if b.index(x) in indexes else 0, b))
            student_checks = check(student_b, table)

            fisher = fisher_kriteria(table, student_checks, SB)
            fisher_check = fisher[0]
            Fp = fisher[1]
            '''print('{} ітерація'.format(counter))
            print('Fp:')
            print(Fp)
            print(Gp*Sbeta*Fp)
            print('\n')'''
            if fisher_check:
                return table, b, check_b, student_b, student_checks, Fp
            else:
                flag = 1
                counter += 1
                x1 = list(map(lambda x: x/(Gp*Sbeta*Fp), x1))
                x2 = list(map(lambda x: x/(Gp*Sbeta*Fp), x2))
                x3 = list(map(lambda x: x/(Gp*Sbeta*Fp), x3))
                xcmax = int(1/3*(x1[1] + x2[1] + x3[1]))
                xcmin = int(1/3*(x1[0] + x2[0] + x3[0]))
                ymax = 200 + xcmax
                ymin = 200 + xcmin
                y = [ymax, ymin]

def normalize_table(table, x1, x2, x3):
    for row in table:
        for j in range(3):
            if j == 0:
                if row[j] == x1[0]:
                    row[j] = -1
                else:
                    row[j] = 1
            elif j == 1:
                if row[j] == x2[0]:
                    row[j] = -1
                else:
                    row[j] = 1
            elif j == 2:
                if row[j] == x3[0]:
                    row[j] = -1
                else:
                    row[j] = 1
        row[3] = row[0]*row[1]
        row[4] = row[0]*row[2]
        row[5] = row[1]*row[2]
        row[6] = row[0]*row[1]*row[2]
    return table

def check(koeficients, table):
    checks = [koeficients[0] for i in range(8)]

    for i in range(len(table)):
        for j in range(7):
            checks[i] += koeficients[j+1]*table[i][j]
    return checks

def cochrane_kriteria(table, N = 8):
    sigma = [table[i][-1] for i in range(N)]
    Gt = 0.5157
    Gp = max(sigma)/sum(sigma)
    if Gp < Gt:
        return 1, Gp
    else:
        return 0, Gp

def student_kriteria(table, N = 8, m = 3):
    sigma = [table[i][-1] for i in range(N)]
    y_mid = [table[i][-2] for i in range(N)]
    SB = sum(sigma)/N
    Sb = SB/(N*m)
    Sbeta = sqrt(Sb)

    beta = [0]*N
    for i in range(N):
        beta[0] += 1/N*(y_mid[i]*1)
    for i in range(1, len(beta)):
        for j in range(N):
            beta[i] += 1/N*(y_mid[j]*table[j][i-1])

    ttab = 2.120

    t = [0]*N
    for i in range(N):
        t[i] = beta[i]/Sbeta

    indexes = []
    for i in t:
        if i > ttab:
            indexes.append(t.index(i))

    return [indexes, SB, Sbeta]

def fisher_kriteria(table, checks, Sb, N = 8, m = 3, d = 1):
    y_mid = [table[i][-2] for i in range(N)]
    Sad = 0
    for i in range(N):
        Sad += m/(N-d)*(checks[i] - y_mid[i])**2
    Ft = 2.7
    Fp = Sad/Sb
    if Fp > Ft:
        return 0, Fp
    else:
        return 1, Fp



x1 = [15, 45]
print('x1min = {0}, x1max = {1}'.format(x1[0], x1[1]))
x2 = [-35, 15]
print('x2min = {0}, x2max = {1}'.format(x2[0], x2[1]))
x3 = [-35, -5]
print('x3min = {0}, x3max = {1}\n'.format(x3[0], x3[1]))


xcmax = int(1/3*(x1[1] + x2[1] + x3[1]))
xcmin = int(1/3*(x1[0] + x2[0] + x3[0]))
print('Xcmin = {0}, Xcmax = {1}\n'.format(xcmin, xcmax))
ymax = 200 + xcmax
ymin = 200 + xcmin
y = [ymin, ymax]
print('ymin = {0}, ymax = {1}\n'.format(ymin, ymax))

sb_counter = 0
for i in range(100):
    research = experiment(x1, x2, x3, y)
    table = research[0]
    koefs_b = research[1]
    check_b = research[2]
    koefs_sb = research[3]
    for j in range(len(koefs_sb)):
        if koefs_sb[j] != 0:
            sb_counter += 1
    check_sb = research[4]
    Fp = research[5]

print("загальна кількість значимих коефіцієнтів:")
print(sb_counter)
print('\n')

print('Таблиця експерименту')
pprint.pprint(table)
print('\n')
print('b:')
print(koefs_b)
print('\n')
print('перевірка b:')
print(check_b)
print('\n')
print('значимі коефіцієнти:')
print(koefs_sb)
print('\n')
print('перевірка sb:')
print(check_sb)
print('\n')
print('Fp:')
print(Fp)
