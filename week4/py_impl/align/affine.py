
from typing import Tuple
INF_NEG = -10**12

def affine_global(a: str, b: str, match=3, mismatch=-3, gap_open=-5, gap_ext=-1) -> Tuple[int, str, str]:
    n, m = len(a), len(b)
    M = [[-10**12]*(m+1) for _ in range(n+1)]
    X = [[-10**12]*(m+1) for _ in range(n+1)]
    Y = [[-10**12]*(m+1) for _ in range(n+1)]
    tbM = [[0]*(m+1) for _ in range(n+1)]
    tbX = [[0]*(m+1) for _ in range(n+1)]
    tbY = [[0]*(m+1) for _ in range(n+1)]
    M[0][0] = 0
    for i in range(1, n+1):
        X[i][0] = gap_open + (i-1)*gap_ext; tbX[i][0] = 1 if i>1 else 0
    for j in range(1, m+1):
        Y[0][j] = gap_open + (j-1)*gap_ext; tbY[0][j] = 1 if j>1 else 0
    def s(x,y): return match if x==y else mismatch
    for i in range(1, n+1):
        ai = a[i-1]
        for j in range(1, m+1):
            bj = b[j-1]
            open_x   = M[i-1][j] + gap_open
            extend_x = X[i-1][j] + gap_ext
            if open_x >= extend_x: X[i][j] = open_x; tbX[i][j] = 0
            else: X[i][j] = extend_x; tbX[i][j] = 1
            open_y   = M[i][j-1] + gap_open
            extend_y = Y[i][j-1] + gap_ext
            if open_y >= extend_y: Y[i][j] = open_y; tbY[i][j] = 0
            else: Y[i][j] = extend_y; tbY[i][j] = 1
            d = M[i-1][j-1] + s(ai,bj)
            fx = X[i-1][j-1] + s(ai,bj)
            fy = Y[i-1][j-1] + s(ai,bj)
            if d >= fx and d >= fy: M[i][j] = d; tbM[i][j] = 0
            elif fx >= fy: M[i][j] = fx; tbM[i][j] = 1
            else: M[i][j] = fy; tbM[i][j] = 2
    end = [(M[n][m],'M'), (X[n][m],'X'), (Y[n][m],'Y')]
    state = max(end, key=lambda x: x[0])[1]
    i, j = n, m
    qa, ta = [], []
    while i>0 or j>0:
        if state=='M':
            src = tbM[i][j]
            qa.append(a[i-1]); ta.append(b[j-1]); i -= 1; j -= 1
            state = 'M' if src==0 else ('X' if src==1 else 'Y')
        elif state=='X':
            qa.append(a[i-1]); ta.append('-'); i -= 1
            state = 'X' if tbX[i+1][j]==1 else 'M'
        else:
            qa.append('-'); ta.append(b[j-1]); j -= 1
            state = 'Y' if tbY[i][j+1]==1 else 'M'
    best = max(end)[0]
    return best, "".join(reversed(qa)), "".join(reversed(ta))
