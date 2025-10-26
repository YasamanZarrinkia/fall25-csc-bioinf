
# Build with:
#   codon build -release -o bin/align_codon codon_impl/align_codon.py
from sys import argv
from typing import List, Tuple

def read_fasta(path: str) -> str:
    seq: List[str] = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line: continue
            if line.startswith(">"): 
                continue
            seq.append(line)
    return "".join(seq).upper()

def score(a: str, b: str, match: int, mismatch: int) -> int:
    return match if a == b else mismatch

def global_alignment(a: str, b: str, match: int, mismatch: int, gap: int) -> Tuple[int, str, str]:
    n = len(a); m = len(b)
    dp: List[List[int]] = [[0]*(m+1) for _ in range(n+1)]
    tb: List[List[int]] = [[0]*(m+1) for _ in range(n+1)]
    for i in range(1, n+1):
        dp[i][0] = dp[i-1][0] + gap; tb[i][0] = 1
    for j in range(1, m+1):
        dp[0][j] = dp[0][j-1] + gap; tb[0][j] = 2
    for i in range(1, n+1):
        ai = a[i-1]
        for j in range(1, m+1):
            bj = b[j-1]
            s_diag = dp[i-1][j-1] + score(ai,bj,match,mismatch)
            s_up   = dp[i-1][j] + gap
            s_left = dp[i][j-1] + gap
            if s_diag >= s_up and s_diag >= s_left:
                dp[i][j] = s_diag; tb[i][j] = 0
            elif s_up >= s_left:
                dp[i][j] = s_up; tb[i][j] = 1
            else:
                dp[i][j] = s_left; tb[i][j] = 2
    i = n; j = m
    qa: List[str] = []; ta: List[str] = []
    while i>0 or j>0:
        t = tb[i][j]
        if t==0:
            qa.append(a[i-1]); ta.append(b[j-1]); i -= 1; j -= 1
        elif t==1:
            qa.append(a[i-1]); ta.append("-"); i -= 1
        else:
            qa.append("-"); ta.append(b[j-1]); j -= 1
    return dp[n][m], "".join(reversed(qa)), "".join(reversed(ta))

def local_alignment(a: str, b: str, match: int, mismatch: int, gap: int) -> Tuple[int, str, str]:
    n = len(a); m = len(b)
    dp: List[List[int]] = [[0]*(m+1) for _ in range(n+1)]
    tb: List[List[int]] = [[3]*(m+1) for _ in range(n+1)]  # 0 diag,1 up,2 left,3 stop
    best = 0; best_i = 0; best_j = 0
    for i in range(1, n+1):
        ch_a = a[i-1]
        for j in range(1, m+1):
            ch_b = b[j-1]
            s_diag = dp[i-1][j-1] + score(ch_a, ch_b, match, mismatch)
            s_up   = dp[i-1][j] + gap
            s_left = dp[i][j-1] + gap
            v = 0; tbv = 3
            if s_diag > v:
                v = s_diag; tbv = 0
            if s_up > v:
                v = s_up; tbv = 1
            if s_left > v:
                v = s_left; tbv = 2
            dp[i][j] = v; tb[i][j] = tbv
            if v > best:
                best = v; best_i = i; best_j = j
    i = best_i; j = best_j
    qa: List[str] = []; ta: List[str] = []
    while i > 0 and j > 0 and dp[i][j] != 0:
        t = tb[i][j]
        if t == 0:
            qa.append(a[i-1]); ta.append(b[j-1]); i -= 1; j -= 1
        elif t == 1:
            qa.append(a[i-1]); ta.append('-'); i -= 1
        elif t == 2:
            qa.append('-'); ta.append(b[j-1]); j -= 1
        else:
            break
    return best, "".join(reversed(qa)), "".join(reversed(ta))
def semi_global_fitting(q: str, t: str, match: int, mismatch: int, gap: int) -> Tuple[int, str, str]:
    n = len(q); m = len(t)
    dp: List[List[int]] = [[0]*(m+1) for _ in range(n+1)]
    tb: List[List[int]] = [[0]*(m+1) for _ in range(n+1)]
    for i in range(1, n+1):
        dp[i][0] = dp[i-1][0] + gap; tb[i][0] = 1
    for j in range(1, m+1):
        dp[0][j] = 0; tb[0][j] = 2
    for i in range(1, n+1):
        qi = q[i-1]
        for j in range(1, m+1):
            tj = t[j-1]
            s_diag = dp[i-1][j-1] + score(qi,tj,match,mismatch)
            s_up   = dp[i-1][j] + gap
            s_left = dp[i][j-1] + gap
            v = s_diag; tbv = 0
            if s_up > v: v = s_up; tbv = 1
            if s_left > v: v = s_left; tbv = 2
            dp[i][j] = v; tb[i][j] = tbv
    # best position at the end of q (row n), any column j
    jbest = 0; best = dp[n][0]
    for j in range(1, m+1):
        if dp[n][j] > best: best = dp[n][j]; jbest = j
    i = n; j = jbest
    qa: List[str] = []; ta: List[str] = []
    while i>0:
        tdir = tb[i][j]
        if tdir==0:
            qa.append(q[i-1]); ta.append(t[j-1]); i -= 1; j -= 1
        elif tdir==1:
            qa.append(q[i-1]); ta.append("-"); i -= 1
        else:
            qa.append("-"); ta.append(t[j-1]); j -= 1
    while j>0:
        ta.append(t[j-1]); qa.append("-"); j -= 1
    return best, "".join(reversed(qa)), "".join(reversed(ta))

def affine_global(a: str, b: str, match: int, mismatch: int, gap_open: int, gap_ext: int) -> Tuple[int, str, str]:
    INF = -10**12
    n = len(a); m = len(b)
    M: List[List[int]] = [[INF]*(m+1) for _ in range(n+1)]
    X: List[List[int]] = [[INF]*(m+1) for _ in range(n+1)]  # gap in b (up)
    Y: List[List[int]] = [[INF]*(m+1) for _ in range(n+1)]  # gap in a (left)
    tbM: List[List[int]] = [[0]*(m+1) for _ in range(n+1)]
    tbX: List[List[int]] = [[0]*(m+1) for _ in range(n+1)]
    tbY: List[List[int]] = [[0]*(m+1) for _ in range(n+1)]
    M[0][0] = 0
    for i in range(1, n+1):
        X[i][0] = gap_open + (i-1)*gap_ext
        tbX[i][0] = 1 if i>1 else 0
    for j in range(1, m+1):
        Y[0][j] = gap_open + (j-1)*gap_ext
        tbY[0][j] = 1 if j>1 else 0
    def s(x: str, y: str) -> int:
        return match if x==y else mismatch
    for i in range(1, n+1):
        ai = a[i-1]
        for j in range(1, m+1):
            bj = b[j-1]
            # X (up)
            open_x   = M[i-1][j] + gap_open
            extend_x = X[i-1][j] + gap_ext
            if open_x >= extend_x:
                X[i][j] = open_x; tbX[i][j] = 0
            else:
                X[i][j] = extend_x; tbX[i][j] = 1
            # Y (left)
            open_y   = M[i][j-1] + gap_open
            extend_y = Y[i][j-1] + gap_ext
            if open_y >= extend_y:
                Y[i][j] = open_y; tbY[i][j] = 0
            else:
                Y[i][j] = extend_y; tbY[i][j] = 1
            # M
            d  = M[i-1][j-1] + s(ai,bj)
            fx = X[i-1][j-1] + s(ai,bj)
            fy = Y[i-1][j-1] + s(ai,bj)
            if d >= fx and d >= fy:
                M[i][j] = d; tbM[i][j] = 0
            elif fx >= fy:
                M[i][j] = fx; tbM[i][j] = 1
            else:
                M[i][j] = fy; tbM[i][j] = 2
    # choose best end state
    endM = M[n][m]; endX = X[n][m]; endY = Y[n][m]
    best = endM; state = 'M'
    if endX > best: best = endX; state = 'X'
    if endY > best: best = endY; state = 'Y'
    i = n; j = m
    qa: List[str] = []; ta: List[str] = []
    while i>0 or j>0:
        if state == 'M':
            src = tbM[i][j]
            qa.append(a[i-1]); ta.append(b[j-1])
            i -= 1; j -= 1
            state = 'M' if src==0 else ('X' if src==1 else 'Y')
        elif state == 'X':
            qa.append(a[i-1]); ta.append('-')
            state = 'X' if tbX[i][j]==1 else 'M'
            i -= 1
        else:
            qa.append('-'); ta.append(b[j-1])
            state = 'Y' if tbY[i][j]==1 else 'M'
            j -= 1
    return best, "".join(reversed(qa)), "".join(reversed(ta))

def main():
    method = "local"
    qpath = ""; tpath = ""
    match = 3; mismatch = -3; gap = -2; gopen = -5; gext = -1
    i = 1
    while i < len(argv):
        a = argv[i]
        if a == "--method": method = argv[i+1]; i += 2
        elif a == "--query": qpath = argv[i+1]; i += 2
        elif a == "--target": tpath = argv[i+1]; i += 2
        elif a == "--match": match = int(argv[i+1]); i += 2
        elif a == "--mismatch": mismatch = int(argv[i+1]); i += 2
        elif a == "--gap": gap = int(argv[i+1]); i += 2
        elif a == "--gap_open": gopen = int(argv[i+1]); i += 2
        elif a == "--gap_ext": gext = int(argv[i+1]); i += 2
        else: i += 1
    q = read_fasta(qpath); t = read_fasta(tpath)
    if method == "global":
        sc, qa, ta = global_alignment(q,t,match,mismatch,gap)
    elif method == "local":
        sc, qa, ta = local_alignment(q,t,match,mismatch,gap)
    elif method == "fitting":
        sc, qa, ta = semi_global_fitting(q,t,match,mismatch,gap)
    else:
        sc, qa, ta = affine_global(q,t,match,mismatch,gopen,gext)
    print(f"score={sc}")
    if len(qa) <= 200:
        print(qa); print(ta)
    else:
        print("(alignment suppressed)")

if __name__ == "__main__":
    main()
