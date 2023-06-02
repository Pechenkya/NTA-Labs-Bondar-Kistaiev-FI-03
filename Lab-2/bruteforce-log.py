def BruteDL(a, b, p):
    ax = 1
    x = 0
    for i in range(0, p):
        if ax == b: return x

        ax = (ax * a) % p
        x += 1

    return float('nan')


import time
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
        x = BruteDL(a, b, p)
        et = time.time()

        print(f"General execution time: {et - st} seconds\n")
        print(f"Solution: {x}")
        print("----------------------------------")
        input("Press enter to continue")