
from typing import Tuple
from .common import score

def global_alignment(a: str, b: str, match=3, mismatch=-3, gap=-2) -> Tuple[int, str, str]:
    n, m = len(a), len(b)
    dp = [[0]*(m+1) for _ in range(n+1)]
    tb = [[0]*(m+1) for _ in range(n+1)]  # 0 diag,1 up,2 left
    for i in range(1, n+1):
        dp[i][0] = dp[i-1][0] + gap; tb[i][0] = 1
    for j in range(1, m+1):
        dp[0][j] = dp[0][j-1] + gap; tb[0][j] = 2
    for i in range(1, n+1):
        ai = a[i-1]
        for j in range(1, m+1):
            bj = b[j-1]
            s_diag = dp[i-1][j-1] + score(ai, bj, match, mismatch)
            s_up   = dp[i-1][j] + gap
            s_left = dp[i][j-1] + gap
            if s_diag >= s_up and s_diag >= s_left:
                dp[i][j] = s_diag; tb[i][j] = 0
            elif s_up >= s_left:
                dp[i][j] = s_up; tb[i][j] = 1
            else:
                dp[i][j] = s_left; tb[i][j] = 2
    i, j = n, m
    qa, ta = [], []
    while i>0 or j>0:
        t = tb[i][j]
        if t == 0:
            qa.append(a[i-1]); ta.append(b[j-1]); i -= 1; j -= 1
        elif t == 1:
            qa.append(a[i-1]); ta.append("-"); i -= 1
        else:
            qa.append("-"); ta.append(b[j-1]); j -= 1
    return dp[n][m], "".join(reversed(qa)), "".join(reversed(ta))
