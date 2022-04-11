#!/usr/bin/env python
# coding: utf-8

import sympy
import numpy as np
from math import modf


def Difference(f_values, x_vals, c=0, n=10):
    del_values = []
    cal_x = True
    if len(x_vals)==0: cal_x = False
        
    for i in range(len(f_values)-1):
        nr = f_values[i+1] - f_values[i]
        dr = 1
        if cal_x:
            dr = x_vals[i+1+c] - x_vals[i]
        del_values.append(round(nr/dr, n))
    return del_values

def DifferenceArray(f_values, x_vals=[], rn=10):    
    Fwd_diff_table = []
    n = len(f_values)
    c = 0
    while len(f_values)!=0:
        if round(sum(f_values)/len(f_values), 5) == f_values[0]:
            Fwd_diff_table.append(f_values)
            Fwd_diff_table += [['-'] * n] * (n - len(Fwd_diff_table))
            break
            
        Fwd_diff_table.append(f_values)
        f_values = Difference(f_values, x_vals, c, rn)
        c+=1
        
    for i in range(len(Fwd_diff_table)):
        if i!=len(Fwd_diff_table):
            Fwd_diff_table[i] += ['-']*(n - len(Fwd_diff_table[i]))
    
    return Fwd_diff_table

def getLength(Fwd_diff_table):
    c = 0
    for i in Fwd_diff_table:
        if (i.count('-')) != len(Fwd_diff_table[0]):
            c+=1
    if c == len(Fwd_diff_table):
        return c
    else: return c+1

def getOperator(Fwd_diff_table, op):    
    del_string = ''
    for i in range(getLength(Fwd_diff_table)):
        if i>1:
            del_string += '|' + str(f'{op}^'+str(i)+' f(x)').center(19)
        elif i==1: del_string += '|'+ f'{op} f(x)'.center(19)
            
    return del_string

def drawOperator(op, num):
    a = {}
    for i in range(1, num+1):
        if i==1:
            a[i] = f'{op} '
        else:
            a[i] = f'{op}^{i} '
    return a

def fac(n):
    if n==0 or n==1:
        return 1
    else: return n*fac(n-1)

def DrawForwardDiffTable(x, f_values, write=False, rn=20):
    Fwd_diff_table = DifferenceArray(f_values)
    del_string = getOperator(Fwd_diff_table, '∆')
    del_op = drawOperator('∆', getLength(Fwd_diff_table))
    
    txt_file = ''
    dashLen = len(f'| {"x".center(4)} |  {"f(x)".center(6)} {del_string}|')
    
    print('-'*dashLen)
    print(f'| {"x".center(4)} |  {"f(x)".center(6)} {del_string}|')
    print('-'*dashLen)
    t = Fwd_diff_table
    
    txt_file += 'Forward Difference Table'.center(150) + '\n\n\n'
    txt_file += '-'*dashLen+'\n'
    txt_file += f'| {"x".center(4)} |  {"f(x)".center(6)} {del_string}|' + '\n'
    txt_file += '-'*dashLen+'\n'

    
    tt = 0 if getLength(Fwd_diff_table)==len(Fwd_diff_table) else 1
    for j in range(len(Fwd_diff_table)-1):
        temp_string = f'| {str(x[j]).center(4)} |'
        i=0
        for k in range(len(Fwd_diff_table[j])):
            if k==0:
                temp_string += f'{t[i][j]}'.center(9) + '|'
            else:
                elt1, elt2 = t[i][j+1] , t[i][j]
                if k==1:
                    if elt1!='-' and elt2!='-':
                        mod = ''
                        if k+1 != len(Fwd_diff_table[j])-tt: mod = ' |'
                        temp_string += f' {del_op[k]}f({x[j]}) = {round(elt1-elt2, rn)}'.ljust(18) + mod
                    else:
                        temp_string += f'\t {str(t[i][j]).ljust(5)}\t|'.center(18)
                else:
                    if elt1!='-' and elt2!='-':
                        mod = ''
                        if k+1 != len(Fwd_diff_table[j])-tt: mod = ' |'
                        temp_string += f' {del_op[k]}f({x[j]}) = {round(elt1-elt2,rn)}'.ljust(18) + mod
                            
                i+=1

        print(temp_string)
        print('-'*dashLen)

        txt_file += temp_string+'\n'
        txt_file += '-'*dashLen+'\n'

    temp_string = f'|{str(x[-1]).center(6)}|{str(t[0][-1]).center(9)}' + '|'

    print(temp_string)
    print('-'*dashLen)

    txt_file += temp_string + ' |'.rjust(dashLen - 18) + '\n'
    txt_file += '-'*dashLen+'\n'

    if write or len(Fwd_diff_table)>5:
        with open("Forward Difference Table.txt", 'w', encoding="utf-8") as f:
            f.write(txt_file)

def DrawBackDiffTable(x, f_values, write=False, rn=10):
    Fwd_diff_table = DifferenceArray(f_values)
    del_string = getOperator(Fwd_diff_table, '∇')
    del_op_inv = drawOperator('∇', getLength(Fwd_diff_table))
    
    txt_file = ''
    dashLen = len(f'| {"x".center(4)} |{"f(x)".center(6)} {del_string}|') + 2
    
    print('-'*dashLen)
    print(f'| {"x".center(4)} |  {"f(x)".center(6)} {del_string}|')
    print('-'*dashLen)
    t = Fwd_diff_table

    txt_file += 'Backward Difference Table'.center(150) + '\n\n'
    txt_file += '-'*dashLen +'\n'
    txt_file += f'| {"x".center(4)} |  {"f(x)".center(6)} {del_string}|' + '\n'
    txt_file += '-'*dashLen +'\n'

    c = 0

    for j in range(len(Fwd_diff_table)):
        temp_string = f'| {str(x[j]).center(4)} |'
        i, l = 0, 0
        tt = 0 if getLength(Fwd_diff_table)==len(Fwd_diff_table) else 1
        try:
            for k in range(len(Fwd_diff_table[j])):
                if k==0:
                    temp_string += f'{t[i][j]}'.center(9) + '|'
                else:
                    elt1, elt2 = t[i][j+1] , t[i][j]

                    if elt1!='-' and elt2!='-':
                        mod = ''
                        if k+1 != len(Fwd_diff_table[j])-tt: mod = ' |'
                        temp_string += f' {del_op_inv[k]}f({x[c:][l+1]}) = {round(elt1-elt2,rn)}'.ljust(18) + mod
            
                    i+=1
                    l+=1

            c+=1        
            print(temp_string)
            print('-'*dashLen)

            txt_file += temp_string
            txt_file += '-'*dashLen +'\n'

            x_temp = x_temp[c:]

        except:pass
    
    temp_string = f'| {str(x[-1]).center(4)} |{str(t[0][-1]).center(9)}' + '|'

    print(temp_string)
    print('-'*dashLen)

    txt_file += temp_string +  ' |'.rjust(dashLen - 18) +'\n'
    txt_file += '-'*dashLen+'\n'

    if write or len(Fwd_diff_table)>6:
        with open("Backward Difference Table.txt", 'w', encoding="utf-8") as f:
            f.write(txt_file)
            
def lgrange_fmlua_UI(x_list, f_values):
    eqns = []
    d = []
    n = []

    print("\t Using Lagrange's Formula for Unequal Interval \t\n\nf(x) = ")
    for i in range(len(x_list)):
        tempEq = []
        denon = 1
        denon_str = ''
        num_str = ''

        for j in range(len(x_list)):
            if i!=j:
                num_str += f'(x- {x_list[j]}) '
                denon *= (x_list[i] - x_list[j])
                denon_str += f'({x_list[i]} - {x_list[j]})'
                tempEq.append([1, -x_list[j]])
        print(f'    {num_str}/({denon_str}) * {f_values[i]}') # f({x_list[i]}) = , denominator: ({denon_str}) 
        d.append(denon)
        n.append(f_values[i])
        eqns.append(tempEq)

    eqns2 = []
    for i in range(len(eqns)):
        p = np.poly1d([1])
        for j in range(len(eqns[i])):
            p *= np.poly1d(eqns[i][j])
        eqns2.append(p)

    print('\n=  ',end=' ')
    p = np.poly1d([0])
    for i in range(len(eqns2)):
        meq = eqns2[i].coeffs
        nr_eqn = " ".join(coef_2_str(meq))
        if i:print('     ',end='')
            
        print(f'({str_to_expr(nr_eqn)}) / {d[i]} * {n[i]}')
        p += meq * n[i]/d[i]

    expr = " ".join(coef_2_str(p.coef[::-1]))
    print('\n=   ' + str(str_to_expr(expr)))

def Newton_Diff_Table(x_values, f_values, bkwd=False):
    Fwd_diff_table = DifferenceArray(f_values)

    exp = input('Which value you want to find: ')
    try:
        value_to_find = float(exp)
    except:
        value_to_find = symbols(exp)

    op = '-'
    a = x_values[0]

    if bkwd:
        a = x_values[-1]
        op = '+'

    h = round(x_values[1] - x_values[0], 5)
    u = check_d_n_r((value_to_find-a)/h)
    n = len(x_values)-1

    cal = 0
    j = 0
    title = "Backward " if bkwd else "Forward "
    print("\nUsing Newton's " + title + " Difference Formula.")
    print(f'\n{"b" if bkwd else "a"}: {a}, h: {h}, u = ({value_to_find} - {a})/{h} : {u}\n')

    var = '-1-j' if bkwd else '0'
    while n+1:
        prod_cal = 1
        c = 0

        if j>0:
            for i in range(j):
                print(f'({u} {op} {c}) ', end=' ')
                prod_cal *=  eval(f'{u}{op}{c}')
                c+=1
        
        if Fwd_diff_table[j][eval(var)] == '-':
            break
            
        if not j: temp = ''
        else: temp = f'/ {j}! * '
            
        print(f'{temp}{Fwd_diff_table[j][eval(var)]} = {check_d_n_r(prod_cal/fac(j)*Fwd_diff_table[j][eval(var)])}')
        cal += prod_cal/fac(j) * Fwd_diff_table[j][eval(var)]
        j+=1; n-=1
        
    cal = check_d_n_r(cal)
    print(f"\n\nThe Value of f({value_to_find}) =  {cal}")
    return cal

def coef_2_str(rx):
    rem = []
    for i in range(len(rx)-1, -1,-1):
        r1 = ''
        if rx[i] >= 0:
            r1 = '+ '
        if i == 0:
            r1 += f'{check_d_n_r(rx[i])}'
        elif i == 1:
            r1 += f'{check_d_n_r(rx[i])}*x'
        elif i != 0:
            r1 +=  f'{check_d_n_r(rx[i])}*x**{i}'
        rem.append(r1)
    return rem

def polyMul(p1, p2):
    return np.poly1d(p1)*np.poly1d(p2)

def str_to_expr(expr):
    return sympy.parsing.sympy_parser.parse_expr(expr)

def check_d_n_r(num, up=5):
    if modf(num)[0]==0:
        return int(num)
    else:
        return round(float(num), up)

def check_decimal(array, n=5):
    ret_arr = []
    for i in range(len(array)):
        ret_arr.append(check_d_n_r(array[i]))
    return ret_arr

def Numerical_Diff(x_values, f_values, show_working=False, up = 2):
    Fwd_diff_table = DifferenceArray(f_values)
    exp = input('Which value you want to find: ')
    value_to_find = check_d_n_r(float(exp))

    a = x_values[0]
    h = round(x_values[1] - x_values[0], 5)
    n = len(x_values)
    x = check_d_n_r(round((value_to_find - a)/h, 5))
    j = 0

    div = ''
    add = '+'
    out = []
    while n:
        temp = ''
        tt = ''
        add2 = ''
        if j>0:
            for i in range(j):
                mul = '*'
                if i: tt += f' {mul} (x - {i})'
                else: tt += f'x'     
                    
        if j: add2 = '+'
        if  j>=1: tt += ' * '
        if j+1 != n and j: div = '/'
        if not j : temp = ''
        else: temp = f'/ {j}!'

        con = f"({sympy.expand('*'.join(tt.split(' * ')[:-1]))}) " if tt else '1'
        con2 = f'/ {fac(j)}' if fac(j)>1 else ''
        
        if Fwd_diff_table[j][0] != '-' and con:
            res = str(sympy.simplify(con)) + f' * {Fwd_diff_table[j][0]} {con2}'
            if res != '0':out.append(res)
        
        if show_working:
            print(add2, tt, end=' ')
            print(f"{Fwd_diff_table[j][0]} {temp}")
            
        j+=1  
        n-=1
        
    if show_working: print('\n = ' + " + ".join(out))

    return out , h, x, value_to_find

def cal_derivative(eq, h, x, at_x, orderOfDerivative=1, up=5):
    prime = "'"*orderOfDerivative
    add_star = f' ** {orderOfDerivative})' if orderOfDerivative>1 else ')'
    t = ''
    
    for i in range(len(eq)):
        add = ''
        if i+1 != len(eq): add = ' + '
        temp_der = sympy.diff(eq[i], 'x', orderOfDerivative)

        if type(temp_der) in [sympy.core.numbers.Zero, sympy.core.numbers.Float]:
            temp_der = check_d_n_r(temp_der)
        
        t += str(temp_der) + add
    
    derivative = sympy.simplify(t)
    f_a = derivative.subs({'x':x})
    
    print(f"\n  Here a = {at_x}, h = {h}, x = {x}")
    print(f"\n ∴ (h{add_star} * f{prime}(a+xh) = ({t})")
    print(f"\n ∴ f{prime}({at_x}) = ({check_d_n_r(f_a)}) / ({h}{add_star}, at x = {x}")
    print(f"\n ∴ f{prime}({at_x}) = {check_d_n_r(f_a/h**orderOfDerivative)}")
    print('--'*30)
