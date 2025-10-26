from typing import Tuple
from .common import score

def local_alignment(a: str, b: str, match=3, mismatch=-3, gap=-2) -> Tuple[int, str, str]:
    n, m = len(a), len(b)
    dp = [[0]*(m+1) for _ in range(n+1)]
    tb = [[3]*(m+1) for _ in range(n+1)]  # 0 diag, 1 up, 2 left, 3 stop
    best, best_i, best_j = 0, 0, 0

    for i in range(1, n+1):
        ch_a = a[i-1]
        for j in range(1, m+1):
            ch_b = b[j-1]
            s_diag = dp[i-1][j-1] + score(ch_a, ch_b, match, mismatch)
            s_up   = dp[i-1][j] + gap           # gap in b (delete b[j-1])
            s_left = dp[i][j-1] + gap           # gap in a (insert b[j-1])

            v = 0
            t = 3  # stop by default when v==0
            if s_diag > v:
                v, t = s_diag, 0
            if s_up > v:
                v, t = s_up, 1
            if s_left > v:
                v, t = s_left, 2

            dp[i][j] = v
            tb[i][j] = t
            if v > best:
                best, best_i, best_j = v, i, j

    # Traceback from best cell until we hit zero
    i, j = best_i, best_j
    qa, ta = [], []
    while i > 0 and j > 0 and dp[i][j] > 0:
        t = tb[i][j]
        if t == 0:      # diag
            qa.append(a[i-1]); ta.append(b[j-1]); i -= 1; j -= 1
        elif t == 1:    # up (gap in b)
            qa.append(a[i-1]); ta.append("-"); i -= 1
        elif t == 2:    # left (gap in a)
            qa.append("-"); ta.append(b[j-1]); j -= 1
        else:
            break

    return best, "".join(reversed(qa)), "".join(reversed(ta))
