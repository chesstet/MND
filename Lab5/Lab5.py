import random as rn
from math import sqrt
from scipy.stats import f, t
import pprint
from sklearn import linear_model

def exp_row(x1_num, x2_num, x3_num, y, m = 3):
    y_gen = [y[0] + (y[1]-y[0])*rn.random() for i in range(m)]
    y_mid = sum(y_gen)/m
    sigma = 0
    for i in range(len(y_gen)):
        sigma += ((y_gen[i] - y_mid)**2)/m
    return [x1_num, x2_num, x3_num, x1_num*x2_num, x1_num*x3_num, x2_num*x3_num, x1_num*x2_num*x3_num, x1_num*x1_num, x2_num*x2_num, x3_num*x3_num, y_gen, y_mid, sigma]

def experiment(x1, x2, x3, y, m = 3, N = 15):
    flag = 1
    counter = 1
    while flag:
        flag = 0

        l = 1.215
        x01 = sum(x1)/2
        x02 = sum(x2)/2
        x03 = sum(x3)/2
        deltax1 = x1[1]-x01
        deltax2 = x2[1]-x02
        deltax3 = x3[1]-x03

        exp1 = exp_row(x1[0], x2[0], x3[0], y)
        exp2 = exp_row(x1[0], x2[0], x3[1], y)
        exp3 = exp_row(x1[0], x2[1], x3[0], y)
        exp4 = exp_row(x1[0], x2[1], x3[1], y)
        exp5 = exp_row(x1[1], x2[0], x3[0], y)
        exp6 = exp_row(x1[1], x2[0], x3[1], y)
        exp7 = exp_row(x1[1], x2[1], x3[0], y)
        exp8 = exp_row(x1[1], x2[1], x3[1], y)
        exp9 = exp_row(-l*deltax1 + x01, x02, x03, y)
        exp10 = exp_row(l*deltax1 + x01, x02, x03, y)
        exp11 = exp_row(x01, -l*deltax2 + x02, x03, y)
        exp12 = exp_row(x01, l*deltax2 + x02, x03, y)
        exp13 = exp_row(x01, x02, -l*deltax3 + x03, y)
        exp14 = exp_row(x01, x02, l*deltax3 + x03, y)
        exp15 = exp_row(x01, x02, x03, y)

        table = [exp1, exp2, exp3, exp4, exp5, exp6, exp7, exp8, exp9, exp10, exp11, exp12, exp13, exp14, exp15]
        norm_table = normalize_table(table, x1, x2, x3)


        cochrane = cochrane_kriteria(table, N)
        cochrane_check = cochrane[0]
        Gp = cochrane[1]

        if not cochrane_check:
            flag = 1
            continue
        else:

            aver_y = list(map(lambda x: x[-2], table))
            x_list = list(map(lambda x: x[0:10], table))
            for i in x_list:
                i.insert(0, 1)

            skm = linear_model.LinearRegression(fit_intercept=False)
            skm.fit(x_list, aver_y)
            b = list(skm.coef_)


            check_b = check(b, table)

            student = student_kriteria(table, N)
            indexes = student[0]
            SB = student[1]
            Sbeta = student[2]
            student_b = list(map(lambda x: x if b.index(x) in indexes else 0, b))
            student_checks = check(student_b, table)

            fisher = fisher_kriteria(table, student_checks, SB)
            fisher_check = fisher[0]
            Fp = fisher[1]
            print('{} ітерація'.format(counter))
            print('Fp:')
            print(Fp)
            print(Gp*Sbeta*Fp)
            print('\n')

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
    l = 1.215
    x01 = sum(x1)/2
    x02 = sum(x2)/2
    x03 = sum(x3)/2

    for i in range(len(table)):
        for j in range(3):
            if i < 8:
                if j == 0:
                    if table[i][j] == x1[0]:
                        table[i][j] = -1
                    else:
                        table[i][j] = 1
                elif j == 1:
                    if table[i][j] == x2[0]:
                        table[i][j] = -1
                    else:
                        table[i][j] = 1
                elif j == 2:
                    if table[i][j] == x3[0]:
                        table[i][j] = -1
                    else:
                        table[i][j] = 1
            else:
                if j == 0:
                    if table[i][j] == -l*(x1[1]-x01) + x01:
                        table[i][j] = -l
                    elif table[i][j] == l*(x1[1]-x01) + x01:
                        table[i][j] = l
                    else:
                        table[i][j] = 0
                elif j == 1:
                    if table[i][j] == -l*(x2[1]-x02) + x02:
                        table[i][j] = -l
                    elif table[i][j] == l*(x2[1]-x02) + x02:
                        table[i][j] = l
                    else:
                        table[i][j] = 0
                elif j == 2:
                    if table[i][j] == -l*(x3[1]-x03) + x03:
                        table[i][j] = -l
                    elif table[i][j] == l*(x3[1]-x03) + x03:
                        table[i][j] = l
                    else:
                        table[i][j] = 0
        table[i][3] = table[i][0]*table[i][1]
        table[i][4] = table[i][0]*table[i][2]
        table[i][5] = table[i][1]*table[i][2]
        table[i][6] = table[i][0]*table[i][1]*table[i][2]
        table[i][7] = table[i][0]*table[i][0]
        table[i][8] = table[i][1]*table[i][1]
        table[i][9] = table[i][2]*table[i][2]
    return table

def check(koeficients, table):
    checks = [koeficients[0] for i in range(15)]

    for i in range(len(table)):
        for j in range(10):
            checks[i] += koeficients[j+1]*table[i][j]
    return checks

def cochrane_kriteria(table, N = 15):
    sigma = [table[i][-1] for i in range(N)]
    Gt = 0.3346
    Gp = max(sigma)/sum(sigma)
    if Gp < Gt:
        return 1, Gp
    else:
        return 0, Gp

def student_kriteria(table, N = 15, m = 3):
    sigma = [table[i][-1] for i in range(N)]
    y_mid = [table[i][-2] for i in range(N)]
    SB = sum(sigma)/N
    Sb = SB/(N*m)
    Sbeta = sqrt(Sb)

    beta = [0]*11
    for i in range(N):
        beta[0] += 1/N*(y_mid[i]*1)
    for i in range(1, len(beta)):
        for j in range(N):
            beta[i] += 1/N*(y_mid[j]*table[j][i-1])

    ttab = 2.042

    t = [0]*11
    for i in range(11):
        t[i] = beta[i]/Sbeta

    indexes = []
    for i in t:
        if i > ttab:
            indexes.append(t.index(i))

    return [indexes, SB, Sbeta]

def fisher_kriteria(table, checks, Sb, N = 15, m = 3, d = 4):
    y_mid = [table[i][-2] for i in range(N)]
    Sad = 0
    for i in range(N):
        Sad += m/(N-d)*(checks[i] - y_mid[i])**2
    Ft = 2.16
    Fp = Sad/Sb
    if Fp > Ft:
        return 0, Fp
    else:
        return 1, Fp



x1 = [-10, 3]
print('x1min = {0}, x1max = {1}'.format(x1[0], x1[1]))
x2 = [-7, 2]
print('x2min = {0}, x2max = {1}'.format(x2[0], x2[1]))
x3 = [-1, 6]
print('x3min = {0}, x3max = {1}\n'.format(x3[0], x3[1]))


xcmax = int(1/3*(x1[1] + x2[1] + x3[1]))
xcmin = int(1/3*(x1[0] + x2[0] + x3[0]))
print('Xcmin = {0}, Xcmax = {1}\n'.format(xcmin, xcmax))
ymax = 200 + xcmax
ymin = 200 + xcmin
y = [ymin, ymax]
print('ymin = {0}, ymax = {1}\n'.format(ymin, ymax))

research = experiment(x1, x2, x3, y)
table = research[0]
xses = list(map(lambda x: x[0:10], table))
yeks = list(map(lambda x: x[10:-1], table))

koefs_b = research[1]
check_b = research[2]
koefs_sb = research[3]
check_sb = research[4]
Fp = research[5]


print('Таблиця експерименту\n')
print("Таблиця іксів")
pprint.pprint(xses)
print("\n")
print("Таблиця ігриків")
pprint.pprint(yeks)
#pprint.pprint(table)
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
