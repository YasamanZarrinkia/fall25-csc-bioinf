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

def upgma(distances):
    D0 = _as_float_matrix(distances)
    _check_sym_nonneg(D0)
    n = len(D0)
    if n == 0:
        raise ValueError("Empty distance matrix")

    nodes = [TreeNode(index=i) for i in range(n)]
    active: List[bool] = [True for _ in range(n)]
    size:   List[int]  = [1 for _ in range(n)]
    height: List[float]= [0.0 for _ in range(n)]
    D: List[List[float]] = [[D0[i][j] for j in range(n)] for i in range(n)]

    merges = 0
    last_i = 0

    while merges < n - 1:
        dist_min = float('inf')
        i_min = -1
        j_min = -1
        for i in range(n):
            if not active[i]:
                continue
            for j in range(i):
                if not active[j]:
                    continue
                dij = D[i][j]
                if dij < dist_min:
                    dist_min = dij
                    i_min = i
                    j_min = j
        if i_min == -1:
            break

        h = dist_min * 0.5
        li = h - height[i_min]
        lj = h - height[j_min]
        nodes[i_min] = TreeNode([nodes[i_min], nodes[j_min]], [float(li), float(lj)])
        height[i_min] = h
        active[j_min] = False

        # update distances to new cluster
        for k in range(n):
            if not active[k] or k == i_min:
                continue
            mean = (
                D[i_min][k] * float(size[i_min]) +
                D[j_min][k] * float(size[j_min])
            ) / float(size[i_min] + size[j_min])
            D[i_min][k] = mean
            D[k][i_min] = mean

        size[i_min] = size[i_min] + size[j_min]
        merges += 1
        last_i = i_min

    # find remaining active as root
    root_idx = last_i
    for i in range(n):
        if active[i]:
            root_idx = i
    return Tree(nodes[root_idx])
