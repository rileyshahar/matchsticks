from sys import argv
from math import comb

def s(m, n):
    if n == 0:
        return 2 ** m
    if m == 0:
        return 2 ** n

    return s(m, n-1) + sum((
        comb(m+1, k) * s(n-1, k)
        for k in range(0, m+1)
    ))

if __name__ == "__main__":
    print(s(int(argv[1]), int(argv[2])))
