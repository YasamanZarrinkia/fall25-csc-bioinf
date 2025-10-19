from typing import List
from tree import Tree, TreeNode

def _as_float_matrix(a: List[List[float]]):
    n = len(a)
    mtx: List[List[float]] = []
    for i in range(n):
        ai = a[i]
        if len(ai) != n:
            raise ValueError("Distance matrix must be square")
        row: List[float] = []
        for j in range(n):
            row.append(float(ai[j]))
        mtx.append(row)
    return mtx

def _check_sym_nonneg(mtx: List[List[float]]):
    n = len(mtx)
    for i in range(n):
        if len(mtx[i]) != n:
            raise ValueError("Distance matrix must be square")
        for j in range(i):
            a = mtx[i][j]
            b = mtx[j][i]
            if a < 0.0 or b < 0.0:
                raise ValueError("Distances must be positive")
            if abs(a - b) > 1e-8:
                raise ValueError("Distance matrix must be symmetric")

def neighbor_joining(distances: List[List[float]]):
    D = _as_float_matrix(distances)
    _check_sym_nonneg(D)
    n = len(D)
    if n < 4:
        raise ValueError("At least 4 nodes are required")

    nodes = [TreeNode(index=i) for i in range(n)]
    active: List[bool] = [True for _ in range(n)]
    remaining = n

    while remaining > 2:
        # Special-case: finish with a 3-way root (Biotite-style)
        if remaining == 3:
            last = [idx for idx in range(n) if active[idx]]
            # deterministic order to make Newick stable
            last.sort()
            i, j, k = last[0], last[1], last[2]
            dij = D[i][j]; dik = D[i][k]; djk = D[j][k]
            li = 0.5 * (dij + dik - djk)
            lj = 0.5 * (dij + djk - dik)
            lk = 0.5 * (dik + djk - dij)
            root = TreeNode([nodes[i], nodes[j], nodes[k]],
                            [float(li), float(lj), float(lk)])
            return Tree(root)

        # r[i] = sum of distances to other active nodes
        r: List[float] = [0.0 for _ in range(n)]
        for i in range(n):
            if not active[i]:
                continue
            s = 0.0
            for k in range(n):
                if not active[k] or k == i:
                    continue
                s += D[i][k]
            r[i] = s

        m = float(remaining)
        qmin = float('inf')
        imin = -1; jmin = -1
        for i in range(n):
            if not active[i]:
                continue
            for j in range(i):
                if not active[j]:
                    continue
                q = (m - 2.0) * D[i][j] - r[i] - r[j]
                if q < qmin:
                    qmin = q; imin = i; jmin = j

        i = imin; j = jmin
        if i == -1:
            break  # should not happen with a valid matrix

        li = 0.5 * (D[i][j] + (r[i] - r[j]) / (m - 2.0))
        lj = D[i][j] - li

        # Merge j into i
        nodes[i] = TreeNode([nodes[i], nodes[j]], [float(li), float(lj)])
        active[j] = False

        # Update distances from merged i to all other active k
        for k in range(n):
            if not active[k] or k == i:
                continue
            dik = 0.5 * (D[i][k] + D[j][k] - D[i][j])
            D[i][k] = dik
            D[k][i] = dik

        remaining -= 1
