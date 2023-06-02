
import numpy as np
import time
from math import exp, sqrt, log, gcd, e
import random
# Use factorization from Lab 1
from factor_module import general_factor, check_prime


# ## Генеруємо факторну базу
lb = 2

def gen_factor_base(n):
    c = 3.38 # recommended constant
    B_lim = int(c * exp(0.5 * sqrt(log(n, lb) * log(log(n, lb), lb))))
    print(f"Factor Base limit: {B_lim}")
    
    S = []
    for a in range(2, B_lim):
        if check_prime(a):
            S.append(a)
    
    return S

# ## Генеруємо систему порівнянь 
C = 20

def is_smooth(val, S):
    D = [0] * len(S)
    for i in range(len(S)):
        while val % S[i] == 0:
            D[i] += 1
            val //= S[i]
        
    if val != 1:
        return [False, []]
    else:
        return [True, D]


def gen_equations(a, n, S, number_of_eq=0):
    st = time.time()
    if number_of_eq == 0:
        number_of_eq = len(S) + C

    A = []
    b = []

    curr_power = 1
    curr_val = a

    while len(b) < number_of_eq:
        [smooth, pows] = is_smooth(curr_val, S)
        if smooth:
            A.append(pows)
            b.append(curr_power)

        curr_val = (curr_val * a) % n
        curr_power += 1

        if curr_val == 1:
            break
    et = time.time()

    print(f"Generated equations: {len(b)}")
    print(f"Generating time: {et - st} seconds\n")
    return (np.array(A), np.array(b))


from multiprocessing import Process, Manager
Subprocesses = 2

def equation_subprocess(a, n, S, start_power, A_shared, b_shared, number_of_eq):
    print("Started subprocess")
    curr_power = start_power
    curr_val = pow(a, start_power, n)

    while len(b_shared) < number_of_eq:
        [smooth, pows] = is_smooth(curr_val, S)
        if smooth:
            A_shared.append(pows)
            b_shared.append(curr_power)

        curr_val = (curr_val * a) % n
        curr_power += 1

        if curr_val == 1:
            break

def gen_equations_parallel(a, n, S, number_of_eq=0):
    st = time.time()
    if number_of_eq == 0:
        number_of_eq = len(S) + C

    with Manager() as manager:
        A_shared = manager.list()
        b_shared = manager.list()

        pow_step = (n-1) // Subprocesses
        pow_start = 1
        p_handles = []
        # Start processes
        for _ in range(Subprocesses):
            p_handles.append(Process(target=equation_subprocess, args=(a, n, S, pow_start, A_shared, b_shared, number_of_eq)))
            p_handles[-1].run()
            p_handles[-1].start()
            pow_start += pow_step

        for i in range(Subprocesses):
            p_handles[i].join()
        et = time.time()

        print(f"Generated equations: {len(b_shared)}")
        print(f"Generating time: {et - st} seconds\n")
        return (np.array(A_shared), np.array(b_shared))

# ## Розв'язування системи лінійних порівнянь за модулем
def solve_modular_eq(A_in, b, mod):
    # Create system of Linear Diophante eq by appending new unknowns
    m = len(A_in[0])
    n = len(A_in)
    A = np.concatenate([A_in, np.transpose([b])], axis=1, dtype='object')

    chosen = []

    additional_shift = {}
    # Gauss elimination
    for j in range(m):
        found = False
        for i in range(n):
            if chosen.count(i) != 0:
                continue

            if gcd(A[i][j], mod) == 1:
                chosen.append(i)
                found = True

                inv = pow(A[i][j], -1, mod)
                A[i] = A[i] * inv % mod
                
                for k in range(n):
                    if k != i and A[k][j] != 0:
                        A[k] = (A[k] - A[k][j]*A[i]) % mod
                
                break
        
        if not found:
            for i in range(n):
                if chosen.count(i) != 0 or A[i][j] == 0:
                    continue

                d = gcd(A[i][j], mod)
                if np.all(A[i] % d == 0):
                    chosen.append(i)
                    A[i] //= d

                    inv = pow(A[i][j], -1, mod)
                    A[i] = A[i] * inv % mod
                    for k in range(n):
                        if k != i and A[k][j] != 0:
                            A[k] = (A[k] - A[k][j]*A[i]) % mod
                    break

    print(A)
    solution = []
    for j in range(m):
        added = False
        for i in range(n):
            if A[i][j] != 0:
                solution.append(A[i][-1])
                added = True
                break
        if not added:
            solution.append(0)
    


    return np.array(solution)

# ## Знаходимо відповідний $\log_\alpha \beta$
def find_index_sequential_comparing(a, beta, n, S, A, b):
    curr_ind = 0
    curr_val = beta

    A_list = A.tolist()
    for _ in range(n-1):
        [smooth, pows] = is_smooth(curr_val, S)
        if smooth and A_list.count(pows) != 0:
            corr_a_ind = b[A_list.index(pows)]
            return (corr_a_ind - curr_ind) % (n - 1)
        
        curr_ind += 1
        curr_val = (curr_val * a) % n

    raise RuntimeError("Can't find index!")

def find_index(a, beta, n, S, S_idxs):
    curr_ind = 0
    curr_val = beta

    for _ in range(n-1):
        [smooth, pows] = is_smooth(curr_val, S)
        if smooth:
            corr_a_ind = np.dot(S_idxs, pows)
            return (corr_a_ind - curr_ind) % (n - 1)
        
        curr_ind += 1
        curr_val = (curr_val * a) % n

    raise RuntimeError("Can't find index!")

# ## Загальний алгоритм розв'язання
def solve_brute(alpha, beta, n):
    Base = gen_factor_base(n)
    A, b = gen_equations(alpha, n, Base)
    res = find_index_sequential_comparing(alpha, beta, n, Base, A, b)
    return res

def solve(alpha, beta, n):
    Base = gen_factor_base(n)
    A, b = gen_equations(alpha, n, Base)
    S_idxs = solve_modular_eq(A, b, n - 1)
    res = find_index(alpha, beta, n, Base, S_idxs)
    return res


import os
clear = lambda: os.system('cls')
clear()

if __name__ == "__main__":
    
    while True:
        clear()
        o = str(input('''Methods:
0 -- slow but reliable (no linear equations)
1 -- fast and furious (originally expected method)
q -- quit
Select method: '''))

        if o == 'q':
            break

        if o != '0' and o != '1':
            continue

        print("Enter parameters: ")
        a = int(input("a = "))
        b = int(input("b = "))
        n = int(input("n = "))

        st = time.time()
        if o == '0':
            res = solve_brute(a, b, n)
        elif o == '1':
            res = solve(a, b, n)
        et = time.time()

        print(f"General execution time: {et - st} seconds\n")
        print(f"Solution: {res}")
        input("Press enter to continue")


