import math
import random
import numpy as np


OPTIMUS_PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]
R = {}

for d in OPTIMUS_PRIMES:
    R[d] = [1]
    while R[d].count(R[d][-1]) < 2:
        R[d].append((R[d][-1] * 2) % d)
    R[d].pop()

# ### Метод пробних ділень
def petod_drobnyx_mylen(num):
    b = bin(num)[:1:-1]
    
    if b[0] == '0':
        return 2
    
    for d in OPTIMUS_PRIMES[1::]:
        sum = 0
        for i in range(len(b)):
            sum += int(b[i]) * R[d][i % len(R[d])]
            sum %= d
        
        if sum == 0:
            return d

    return 1

# ### $\rho$-метод Поларда для знаходження дільників
def rho_f(x, mod):
    return (x**2 + 1) % mod

def rho_pollard(num, x_0 = 2, f = rho_f):
    X = [f(x_0, num)]
    Y = [f(X[0], num)]

    i = 0
    while X.count(X[i]) != 2:
        X.append(f(X[i], num))
        Y.append(f(f(Y[i], num), num))
        i += 1

        # print(f"{X[i], Y[i]}")
        if X[i] == Y[i]:
            return 1

        d = math.gcd(num, (Y[i] - X[i]) % num)
        if d != 1:
            return d

    return 1

def jacobi(a, n):
    a = int(a)
    
    if math.gcd(a, n) != 1:
        return 0

    if a == 1:
        return 1

    if a > n:
        return jacobi(a % n, n)

    if a % 2 == 0:
        if n % 8 == 1 or n % 8 == 7:
            return jacobi(a / 2, n)
        else:
            return (-1) * jacobi(a / 2, n)

    if n % 4 == 1 or a % 4 == 1:
        return jacobi(n, a)
    else:
        return (-1) * jacobi(n, a)
    
# ### Ймовірнісний алгоритм Міллера-Рабіна та загальний алгоритм для знаходження простих чисел
def miller_rabin(num, base):
    i = 1
    while (num - 1) % (2 ** i) == 0:
        i += 1

    k = i - 1
    d = (num - 1) // (2 ** k)

    a_d = pow(base, d, num)

    if a_d == 1:
        return True
    
    a_d2i = a_d
    for j in range(k):
        if a_d2i == (num - 1):
            return True
        
        a_d2i = (a_d2i ** 2) % num

    return False


def check_prime(num, error_prob = 0.01):
    if OPTIMUS_PRIMES.count(num) != 0:
        return True

    t = int(math.ceil(math.log(1 / error_prob, 4)))
    s = 0
    for _ in range(t):
        a = random.randrange(3, num + 1)
        s += int(miller_rabin(num, a))

    return s > (t / 2)

# ### Швидкий метод Гауса для знаходження розв'язків СЛР над $GF_2$
def add_columns(A, fr, to):
    if fr == to:
        return      # Smort optimi3aцion
    
    for i in range(len(A)):
        A[i][to] ^= A[i][fr]


# A - square matrix n x m
# returns rows with linear independed vectors
def gaussian_elim_gf2(A):
    m = len(A[0])
    n = len(A) - m
    marked = []                 # marked rows
    for j in range(m):          # cols
        for i in range(n):      # rows
            if A[i][j] == 1:
                marked.append(i)
                # eliminate 1s in this row
                for k in range(m):  
                    if A[i][k] == 1:
                        add_columns(A, j, k)
                break
    return marked

def get_ker_basis(A):
    n = len(A)
    m = len(A[0])
    # build matrix for Gauss elimination
    I = np.identity(m, dtype=int)
    B = np.vstack((A, I)) 
    _marked = gaussian_elim_gf2(B)
    B = np.transpose(B)

    ker_basis = []
    cmp_list = [0 for _ in range(n)]
    for i in range(m):
        if np.all(B[i][0:n:] == cmp_list):
            ker_basis.append(B[i][n::].tolist())

    return ker_basis

# ### Функції знаходження розв'язків квадратного порівнняння і квадратного кореня за модулем
def sqeq_modp(a, b, c, p): # сквек
    k = (b**2 - 4*a*c) % p

    if jacobi(k, p) != 1:
        print(f"error k = {k}, p = {p}")

    y = sqrt_modp(k, p)

    x = [0, 0]
    x[0] = (pow(2*a, -1, p) * (y[0] - b)) % p
    x[1] = (pow(2*a, -1, p) * (y[1] - b)) % p 

    return x


def sqrt_modp(a, p):
    if jacobi(a, p) != 1:
        print(f"error a = {a}, p = {p}")

    if p % 4 == 3:
        # print("4k + 3")
        sq_a = pow(a, (p + 1) // 4, p)
        return [sq_a, p - sq_a]
    
    if p % 8 == 5:
        # print("8k + 5")
        k = (p - 5) // 8
        if pow(a, 2*k + 1, p) == 1:
            sq_a = pow(a, k + 1, p)
        else:
            sq_a = (pow(a, k + 1, p) * pow(2, 2*k + 1, p)) % p

        return[sq_a, p - sq_a]
    
    if p % 8 == 1:
        # print("8k + 1")
        b = 2
        while jacobi(b, p) != -1:
            b = random.randrange(3, p - 1)

        t_a = (p - 1) // 2
        t_b = 0

        while t_a % 2 == 0:
            if (pow(a, t_a, p) * pow(b, t_b, p)) % p  == p - 1:
                t_b += (p - 1) // 2

            t_a = t_a // 2
            t_b = t_b // 2

        if (pow(a, t_a, p) * pow(b, t_b, p)) % p  == p - 1:
                t_b += (p - 1) // 2

        sq_a = (pow(a, (t_a + 1) // 2, p) * pow(b, t_b // 2, p)) % p
        return[sq_a, p - sq_a]

# ### Метод Померанця (Квадратичного Сита)
# 
# Обмеження підібрані під наш варіант: \
# Значення простих чисел: <6000 \
# Інтервал просіювання [-2000, 2000] \
# Максимальний набір гладких чисел: 200
# 
# Такі значення дозволяють знайти найменшу достатню множину розв'язків СЛР (перевірено експериментальним шляхом :)).

Prime_Base_Bound = 10
Sieving_Interval = 100
Max_Smooth = 200


def QS(num, B = Prime_Base_Bound, M = Sieving_Interval, MAX_SMOOTH = Max_Smooth):
    m = math.isqrt(num)
    def q(x):
        return ((x + m)**2 - num) % num
    
    # Prime base formation
    P = []
    for a in range(3, B):
        if check_prime(a):
            if jacobi(num, a) == 1:
                P.append(a)
    logP = [math.log2(p) for p in P]

    # Calculate values in range (optimization)
    Q = []
    for x in range(-M, M):
        q_x = q(x)
        Q.append([x, q_x, math.log2(q_x)])

    # Sieving
    for i in range(len(P)):
        x_1, x_2 = sqeq_modp(1, -2 * m, m**2 - num, P[i])
        for k in range(-(M // P[i]), (M // P[i])):
            pos_x1 = x_1 + k * P[i] + M
            pos_x2 = x_2 + k * P[i] + M
            Q[pos_x1][2] -= logP[i]
            Q[pos_x2][2] -= logP[i]
    Q = sorted(Q, key=lambda a: a[2])

    # Function to check if n is P-smooth
    def smoothness(n):
        powers = [0]*len(P)
        for i in range(len(P)):
            while n % P[i] == 0:
                powers[i] += 1
                n //= P[i]
        if n != 1:
            return []
        return powers
    
    # Find necessary amount of P-smooth numbers
    # Factor this numbers
    Smooth = []
    i = 0
    for [_, qx, _] in Q:
        if i >= MAX_SMOOTH:
            break
        pows = smoothness(qx)
        if len(pows) != 0:
            Smooth.append((qx, pows))
            i += 1

    # Form binary matrix and get basis of it's kernel
    A = [[i % 2 for i in v[1]] for v in Smooth]
    A = np.transpose(A)
    ker_basis = get_ker_basis(A)
    print(f"dim(ker) = {len(ker_basis)}")

    # Utility (add powers to Y representation)
    def add_powers(to, pows):
        for i in range(len(pows)):
            to[i] += pows[i]
    
    # Random-based generator of vectors from ker
    def generate_rand_sol(kernel_vec):
        res = [0] * len(P)
        for v in kernel_vec:
            if random.randint(0, 1) == 1:
                for i in range(len(v)):
                    res[i] += v[i]
        return res

    # Search for divisors
    for _ in range(2**len(ker_basis)):
        sol = generate_rand_sol(ker_basis)

        X = 1
        Y_pows = [0]*len(P)
        for j in range(len(sol)):
            if sol[j] == 1:
                X = (X * Smooth[j][0]) % num
                add_powers(Y_pows, Smooth[j][1])
        Y = 1
        for i in range(len(Y_pows)):
            Y = Y * (P[i]**(Y_pows[i] // 2)) % num
        
        if X == Y:
            continue

        a = math.gcd((X - Y) % num, num)
        if a != 1:
            return a
        
        b = math.gcd((X + Y) % num, num)
        if b != 1:
            return b

    return 1
            

# ### Загальний алгоритм факторизації
def general_factor(num):
    factors = []
    
    # Check if num is prime
    if check_prime(num):
        return [num]
    
    # Check small divs
    while True:
        d = petod_drobnyx_mylen(num)
        if d == 1:
            break
        factors.append(d)
        num //= d

    # rho-Pollard
    d = rho_pollard(num)
    if d != 1:
        factors.append(d)
        num //= d

    # QS till the end
    while True:
        if check_prime(num): 
            factors.append(num)
            return factors

        d = QS(num, B=math.sqrt(num), S = 20, Max_Smooth=math.sqrt(num))
        if d == 1:
            factors.append(num)
            print(f'''Cannot factor {num} for now... X(\n
            Try to change basic parameters''')
            return factors
        num //= d
        factors.append(d)