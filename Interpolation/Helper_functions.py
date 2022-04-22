import sympy
from time import sleep
from math import *
import inspect

def make_table_xy(x_vals, y_vals, f):
    res = []
    expr = sympy.simplify(f)
    for i in range(len(x_vals)):
        temp_res = []
        for j in range(len(y_vals)):
            temp_res.append(expr.subs({'x': x_vals[i], 'y':y_vals[j]}))
        res.append(temp_res)
    return res

def Draw_table(x_vals, y_vals):
    n = max([len(str(i)) for i in x_vals]) + 5
    m = max([len(str(i)) for i in y_vals]) + 5
    
    if m>n: n = m
    x_t = '|' + 'x'.center(n)
    y_t = '|' + 'y'.center(n)
    for i in range(len(x_vals)):
        x_t += '|' + str(x_vals[i]).center(n) 
        
    for j in range(len(y_vals)):
        y_t += '|' + str(y_vals[j]).center(n) 

    x_t += '|'
    y_t += '|'

    if len(x_vals)==len(y_vals):
        print()
        print('-'*len(x_t))
        print(x_t)
        print('-'*len(x_t))
        print(y_t)
        print('-'*len(x_t))
        print()
    else:
        print('[-] Check the Length of X-coords & Y-coords')
    
def create_xy_table(x_vals, y_vals, arr, n=10): 
    temp = "|".join([str(i).center(n) for i in x_vals]) + '  |'
    tt = '|' + 'y \ x'.center(n) + '| ' + temp
    print('\n' + ' ' + '-'*(len(tt)-2))
    print(tt)
    print(' ' + '-'*(len(tt)-2))

    for i in range(len(y_vals)):
        print('|' + str(y_vals[i]).center(n) + '|',end='')
        for j in range(len(x_vals)):
            print(str(arr[j][i]).center(n+1) + '|', end='')
        print('\n' + '-'*len(tt))
    print('\n')
    
def sub_vals(f, *args):
    if len(args)>1:
        return f.subs({'x':args[0], 'y':args[1]})
    else:
        return f.subs({'x':args[0]})
    
    
def cvt_eqn(list, factor, alpha='x'):
    res = []
    for i in range(len(list)):
        add = ''
        mul = ''
        t = f'{list[i]} '
        pwr = '**' + str(len(list)-i-1)
        add_fac = f'({alpha} - {factor})'
        
        if len(list) != i+1: add=' + '
        if len(list)!=i+1 and list[i]!=1: mul = ' * '
        if len(list)-i-1 in [1, 0]: pwr = ''
        if len(list) == i+1: add_fac = ''
        if list[i] == 1: t = ''
        
        a = f'{t}{mul}{add_fac}{pwr}'
        res.append(a)

    return " + ".join(res)


def getFactors(x_vals, l, alpha='x'):
    t  = ''
    for k in range(l):
        mul = ''
        if k+1 != l: mul = ' * '
        c = x_vals[k]
        if c<0: c = f'({c})'
        t += f'({alpha} - {c}){mul}'
    return t


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


def str_to_expr(expr):
    return sympy.parsing.sympy_parser.parse_expr(expr)


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


def getFactors(x_vals, l, alpha='x'):
    t  = ''
    for k in range(l):
        mul = ''
        if k+1 != l: mul = ' * '
        c = x_vals[k]
        if c<0: c = f'({c})'
        t += f'({alpha} - {c}){mul}'
    return t


def synDiv(coef,x_0, printTable = False):
    temp = 0
    out = []
    temp_list = []
    
    for i in coef:
        temp_list.append(round(temp * x_0, 6))
        temp = temp * x_0 + i
        out.append(round(temp, 6))
        
    if out[-1] == 0:
        print(f'{x_0} is the exact root of the given equation')
    
    if printTable:
        print('\n')
        print(f'\t |   ', "\t".join([str(i).rjust(8) for i in coef]))
        print("\t |   ")
        print(f" {str(x_0).ljust(8)}|   ","\t".join([str(i).rjust(8) for i in temp_list]))
        print('----'*25)
        print('\t |   ',"\t".join([str(i).rjust(8) for i in out]))
        print('')

    return out


def getIntervals(x, gap):
    intervals = []
    for i in range(0, len(x), 1):
        if len(x[i:i+gap+1])>gap:
            intervals.append(x[i:i+gap+1])

    return intervals


def get_f_x_y(x_vals, y_vals):
    m, n = len(x_vals), len(y_vals)
    arr = []
    for i in range(m):
        tem_arr = []
        for j in range(n):
            tem_arr.append(float(input(f"value at {x_vals[i], y_vals[j]} : ")))
        tem_arr = check_decimal(tem_arr)
        arr.append(tem_arr)
    return arr


def find_values(arr, x_vals, y_vals, x_find, y_find):
    if x_find in x_vals and y_find in y_vals:
        return arr[x_vals.index(x_find)][y_vals.index(y_find)]
    
    print(f"No value was found at ({x_find}, {y_find})")


def value_at(list, value, i=0):
    return list[list.index(value)+i]


def get_x_y_at(x_vals, y_vals, arr):
    a = {}
    for i in range(len(arr)):
        tset = {}
        for j in range(len(arr[i])):
            tset[y_vals[j]] = arr[i][j]
        a[x_vals[i]] = tset

    return a


def getHeight(x):
    out = []
    for i in range(len(x)-1):
        out.append(round(x[i+1]-x[i], 15))
    if sum(out)/len(out) == out[0]:
        return out[0]
    
    print("The intervals are not equally divided")
    
    
def print_details(x_vals, y_vals, f):
    arr = make_table_xy(x_vals, y_vals, f)
    print(f'f(x, y) = {f}')
    create_xy_table(x_vals, y_vals, arr)
    print('\n')
    
def DrawMat(L, n=8):
    for i in range(len(L)):
        if i==0:
            print(f'{"[" + "".join([str(i).center(n) for i in L[0]]) + "]"}')
        else:
            print(f'          {" ".center(n)} {"[" + "".join([str(i).center(n) for i in L[i]]) + "]"}')
    print('')