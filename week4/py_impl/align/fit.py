
from typing import Tuple
from .common import score

def semi_global_fitting(q: str, t: str, match=3, mismatch=-3, gap=-2) -> Tuple[int, str, str]:
    n, m = len(q), len(t)
    dp = [[0]*(m+1) for _ in range(n+1)]
    tb = [[0]*(m+1) for _ in range(n+1)]
    for i in range(1, n+1):
        dp[i][0] = dp[i-1][0] + gap; tb[i][0] = 1
    for j in range(1, m+1):
        dp[0][j] = 0; tb[0][j] = 2
    for i in range(1, n+1):
        qi = q[i-1]
        for j in range(1, m+1):
            tj = t[j-1]
            s_diag = dp[i-1][j-1] + score(qi, tj, match, mismatch)
            s_up   = dp[i-1][j] + gap
            s_left = dp[i][j-1] + gap
            v = s_diag; tbv = 0
            if s_up > v: v = s_up; tbv = 1
            if s_left > v: v = s_left; tbv = 2
            dp[i][j] = v; tb[i][j] = tbv
    jbest = max(range(m+1), key=lambda j: dp[n][j])
    best = dp[n][jbest]
    i, j = n, jbest
    qa, ta = [], []
    while i>0:
        tdir = tb[i][j]
        if tdir == 0:
            qa.append(q[i-1]); ta.append(t[j-1]); i -= 1; j -= 1
        elif tdir == 1:
            qa.append(q[i-1]); ta.append("-"); i -= 1
        else:
            qa.append("-"); ta.append(t[j-1]); j -= 1
    while j>0:
        ta.append(t[j-1]); qa.append("-"); j -= 1
    return best, "".join(reversed(qa)), "".join(reversed(ta))
