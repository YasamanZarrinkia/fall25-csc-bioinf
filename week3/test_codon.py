# Codon test entry: run ONLY the 3 required tests against the Codon port
import upgma as _upgma
import neighbor_joining as _nj
import tree as _tree
from typing import List

Tree = _tree.Tree
TreeNode = _tree.TreeNode
TreeError = _tree.TreeError
InvalidFileError = _tree.InvalidFileError
upgma = _upgma.upgma
neighbor_joining = _nj.neighbor_joining

def test(f): return f  # identity

# --- pure-Codon loader for distances.txt (returns List[List[float]]) ---
def _load_txt_matrix_float(path: str) -> List[List[float]]:
    rows: List[List[float]] = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            row: List[float] = []
            for p in parts:
                row.append(float(p))
            rows.append(row)
    return rows

def distances_fixture() -> List[List[float]]:
    return _load_txt_matrix_float("testdata/distances.txt")

def upgma_newick_fixture() -> str:
    with open("testdata/newick_upgma.txt", "r") as f:
        return f.read().strip()

# ------------------ tests ------------------

@test
def test_upgma():
    """
    Compare the results of upgma() with the provided Newick (DendroUPGMA).
    """
    dist = distances_fixture()         # List[List[float]]
    tree = upgma(dist)
    ref_newick = upgma_newick_fixture()
    ref_tree = Tree.from_newick(ref_newick)

    n = len(tree)
    for i in range(n):
        for j in range(n):
            d1 = tree.get_distance(i, j)
            d2 = ref_tree.get_distance(i, j)
            assert abs(d1 - d2) <= 1e-3
            assert tree.get_distance(i, j, True) == ref_tree.get_distance(i, j, True)

@test
def test_neighbor_joining():
    # Classic 6x6 NJ example (expects a trifurcating root)
    dist: List[List[float]] = [
        [0.0, 5.0, 4.0, 7.0, 6.0, 8.0],
        [5.0, 0.0, 7.0,10.0, 9.0,11.0],
        [4.0, 7.0, 0.0, 7.0, 6.0, 8.0],
        [7.0,10.0, 7.0, 0.0, 5.0, 9.0],
        [6.0, 9.0, 6.0, 5.0, 0.0, 8.0],
        [8.0,11.0, 8.0, 9.0, 8.0, 0.0],
    ]
    ref_tree = Tree(
        TreeNode(
            [
                TreeNode(
                    [
                        TreeNode([TreeNode(index=0), TreeNode(index=1)], [1.0, 4.0]),
                        TreeNode(index=2),
                    ],
                    [1.0, 2.0],
                ),
                TreeNode([TreeNode(index=3), TreeNode(index=4)], [3.0, 2.0]),
                TreeNode(index=5),
            ],
            [1.0, 1.0, 5.0],
        )
    )
    test_tree = neighbor_joining(dist)
    assert test_tree == ref_tree

@test
def test_distances():
    """
    UPGMA tree must be ultrametric; and specific topological distances should match.
    """
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
    print(f"{'codon':<11}{ms}ms")
