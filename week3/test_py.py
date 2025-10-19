# Python test entry: run ONLY the 3 required tests against Biotite
import numpy as np
import biotite.sequence.phylo as phylo
from biotite.file import InvalidFileError as InvalidFileError

Tree = phylo.Tree
TreeNode = phylo.TreeNode
TreeError = phylo.TreeError
upgma = phylo.upgma
neighbor_joining = phylo.neighbor_joining

def test(f): return f  # identity

# --- pytest.approx replacement (local, minimal) ---
def pytest_approx(value, **kwargs):
    tol = kwargs.get("abs", 1e-12)
    class _Approx:
        def __init__(self, val, tol): self.val, self.tol = val, tol
        def __eq__(self, other):
            import builtins
            return builtins.abs(other - self.val) <= self.tol
    return _Approx(value, tol)

# --- fixtures (same files as Codon side) ---
def distances_fixture():
    return np.loadtxt("testdata/distances.txt", dtype=float)

def upgma_newick_fixture():
    with open("testdata/newick_upgma.txt", "r") as f:
        return f.read().strip()

# ------------------ tests ------------------

@test
def test_upgma():
    dist = distances_fixture()
    tree = upgma(dist)
    ref_newick = upgma_newick_fixture()
    ref_tree = Tree.from_newick(ref_newick)
    n = len(tree)
    for i in range(n):
        for j in range(n):
            assert tree.get_distance(i, j) == pytest_approx(
                ref_tree.get_distance(i, j), abs=1e-3
            )
            assert tree.get_distance(i, j, True) == ref_tree.get_distance(i, j, True)

@test
def test_neighbor_joining():
    dist = np.array([
        [ 0,  5,  4,  7,  6,  8],
        [ 5,  0,  7, 10,  9, 11],
        [ 4,  7,  0,  7,  6,  8],
        [ 7, 10,  7,  0,  5,  9],
        [ 6,  9,  6,  5,  0,  8],
        [ 8, 11,  8,  9,  8,  0],
    ], dtype=float)
    ref_tree = Tree(
        TreeNode(
            [
                TreeNode(
                    [
                        TreeNode([TreeNode(index=0), TreeNode(index=1)], [1, 4]),
                        TreeNode(index=2),
                    ],
                    [1, 2],
                ),
                TreeNode([TreeNode(index=3), TreeNode(index=4)], [3, 2]),
                TreeNode(index=5),
            ],
            [1, 1, 5],
        )
    )
    test_tree = neighbor_joining(dist)
    assert test_tree == ref_tree

@test
def test_distances():
    dist = distances_fixture()
    tree = upgma(dist)
    root = tree.root
    base = root.distance_to(tree.leaves[0])
    for leaf in tree.leaves:
        assert leaf.distance_to(root) == base
    assert tree.get_distance(0, 19, True) == 9
    assert tree.get_distance(4, 2, True) == 10

# --------------- run & timing ---------------
def _now_ms():
    import time
    return time.perf_counter_ns() // 1_000_000

def _run_all():
    test_upgma()
    test_neighbor_joining()
    test_distances()

def _now_ns():
    import time
    return time.perf_counter_ns()


def _run_all_tests_and_time():
    t0 = _now_ns()
    _run_all()
    t1 = _now_ns()
    # ceil(ns / 1e6) without floats
    return (t1 - t0 + 999_999) // 1_000_000

if __name__ == "__main__":
    ms = _run_all_tests_and_time()
    print("Language    Runtime")
    print("-------------------")
    print(f"{'python':<11}{ms}ms")
