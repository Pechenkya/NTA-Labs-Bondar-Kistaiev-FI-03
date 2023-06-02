# %%
import math
import time
from factor_module import check_prime

# %%
# Solves system of congruences x = a[i] mod n[i] by Chinese Remainder Theorem
def CRT(a, n):
    n_prod = math.prod(n)
    N = [n_prod // n_i for n_i in n]
    M = [pow(n_prod // n_i, -1, n_i) for n_i in n]

    return sum([a[i]*M[i]*N[i] for i in range(0, len(a))]) % n_prod


# %%
# Generates all primes from 2 to n
def PrimesRange(n):
    P = [2]
    for a in range(3, n):
        if check_prime(a):
            P.append(a)
    
    return P

def NextPrime(n):
    p = n + 1
    while not check_prime(p):
        p += 1

    return p

# Factors n 
def Factor(n):
    F = {}
    p_i = 2

    while n != 1:
        k = 0
        while n % (p_i ** (k+1)) == 0:
            k += 1

        if k > 0:
            F[p_i] = k
            n = n // (p_i**k)
        
        p_i = NextPrime(p_i)

    if n != 1:
        F[n] = 1

    return F

# %%
# Computes discrete log_a(b) modp by Silver–Pohlig–Hellman algorithm
def SPH(a, b, p):
    n = p - 1
    F = Factor(n)
    r_table = {}

    for p_i in F.keys():
        for j in range(0, p_i):
            t = pow(a, (n * j) // p_i, p)
            r_table[p_i, t] = j
    
    Y = []
    P = [p_i**l_i for p_i, l_i in F.items()]

    for p_i in F.keys():
        a_x_modn = 1
        x = 0
        for j in range(0, F[p_i]):
            t = pow(b * pow(a_x_modn, -1, p), n // (p_i ** (j+1)), p)
            x_j = r_table[p_i, t]
            x += (x_j * (p_i ** j)) % (p_i ** F[p_i])
            a_x_modn = (a_x_modn * pow(a, x_j * (p_i ** j), p)) % p
        
        Y.append(x)
    
    return CRT(Y, P)

# %%
def BruteDL(a, b, p):
    ax = 1
    x = 0
    for i in range(0, p):
        if ax == b: return x

        ax = (ax * a) % p
        x += 1

    return float('nan')

import os
clear = lambda: os.system('clear')

if __name__ == "__main__":
    
    while True:
        clear()
        
        print("----------------------------------")
        print("Enter parameters: ")
        a = int(input("a = "))
        b = int(input("b = "))
        p = int(input("p = "))

        st = time.time()
        x = SPH(a, b, p)
        et = time.time()

        print(f"General execution time: {et - st} seconds\n")
        print(f"Solution: {x}")
        print("----------------------------------")
        input("Press enter to continue")


