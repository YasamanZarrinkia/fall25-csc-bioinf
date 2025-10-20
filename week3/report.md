# Week 3 — Porting Biotite’s `phylo` (UPGMA & Neighbor-Joining) to Codon

> **Acknowledgment**: This report was written with assistance from AI (ChatGPT).

## Goal

Port just enough of Biotite’s `phylo` to **Codon** so these 3 tests pass and are timed:

- `test_upgma`
- `test_neighbor_joining`
- `test_distances`

Each runner prints:

```
Language    Runtime
-------------------
python      <X>ms
codon       <Y>ms
```

**Runtime** = *time needed to run all three tests*, measured **inside** each test program (tests-only; excludes interpreter startup and Codon compile time).

---

## Final layout

```
week3/
  tree.py
  upgma.py
  neighbor_joining.py
  test_codon.py      # Codon tests (the 3 tests only)
  test_py.py         # Python/Biotite reference (same 3 tests)
  evaluate.sh        # runs both and prints one table
  testdata/
    distances.txt
    newick_upgma.txt
.github/
  workflows/
    week3.yml        # CI
```

> The two data files used by the tests go in `week3/testdata/`.

---


## What was ported (Codon)

### 1) `tree.py`
Minimal tree support to match Biotite semantics used by the tests:

- `TreeNode`:  
  `index`, `children`, `parent`, `distance` (root reports **0.0**),  
  `lowest_common_ancestor()`, `distance_to()`, `get_leaves()`, `get_indices()`, `copy()`
- `Tree`: wraps root, fixes leaf order 0..n-1, `get_distance(i,j, topological=False)`, `copy()`, `to_newick()`, `from_newick()`.
- Minimal **Newick parser** sufficient for tests:
  - integer leaf labels 
  - nested groups
  - optional branch lengths after nodes and after groups, whitespace tolerant
  - ignores intermediate labels
- **Equality**: order-insensitive, structural comparison:
  - topological distances must match exactly
  - metric distances must match within small epsilon (e.g., `1e-9`)

### 2) `upgma.py`
Classic UPGMA:

- Active clusters with `size[]` and `height[]`
- Merge the closest pair; new node branch lengths  
  `li = h - height[i]`, `lj = h - height[j]` with `h = dij / 2`
- Update distances to the merged cluster via size-weighted mean

### 3) `neighbor_joining.py`
Standard NJ:

- Active set; matrix `D`
- `r[i] = sum_k D[i,k]` over active
- Minimize `Q = (m-2) D[i,j] − r[i] − r[j]`
- Merge `i,j`; update `dik = 0.5*(D[i,k]+D[j,k]-D[i,j])`
- **Special case when 3 remain** → create a **trifurcating root** (Biotite parity):  
  `li = 0.5*(dij + dik − djk)`, `lj = 0.5*(dij + djk − dik)`, `lk = 0.5*(dik + djk − dij)`  
  (Sort the final 3 indices for deterministic output.)

---

## Tests

### `test_py.py` (Python/Biotite)
- Loads `testdata/distances.txt` with NumPy and `newick_upgma.txt` with `open()`.
- `test_upgma()` compares all pairwise **metric** distances (±1e-3) and **topological** distances to a Newick reference.
- `test_neighbor_joining()` compares to the classic 6×6 NJ reference tree (trifurcating root).
- `test_distances()` checks UPGMA ultrametricity and two sample topological distances.

### `test_codon.py` (Codon)
- Pure Codon line-by-line loader for `distances.txt` to avoid Codon NumPy I/O issues.
- Same 3 tests and assertions as Python.

### Timing (tests-only)
In both files:

- Use `time.perf_counter_ns()` around `_run_all()` (the three tests) and **ceil** ns→ms:
  ```python
  def _run_all_tests_and_time():
      t0 = _now_ns()
      _run_all()
      t1 = _now_ns()
      return (t1 - t0 + 999_999) // 1_000_000  # ceil to ms
  ```
- Print exactly:
  ```
  Language    Runtime
  -------------------
  python      <ms>ms   # in test_py.py
  codon       <ms>ms   # in test_codon.py
  ```


---

## Runner (`evaluate.sh`)

A small wrapper that just executes the two programs and relays their internally measured runtimes:

```bash
#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

PY_OUT="$(python3 test_py.py)"
CD_OUT="$(codon run -release test_codon.py)"

echo "Language    Runtime"
echo "-------------------"
echo "$PY_OUT" | awk '/^python[[:space:]]/{print;exit}'
echo "$CD_OUT" | awk '/^codon[[:space:]]/{print;exit}'
```

---

## CI Workflow (`.github/workflows/week3.yml`)

Installs Python & Biotite, installs Codon + seq plugin, runs `week3/evaluate.sh`, and posts the results in the job summary.

---

## Codon gotchas (and fixes)

1) **No `frozenset`; sets not hashable**  
   Biotite uses `frozenset` in hashing. Add a stable hash for `set`:
   ```python
   @extend
   class set:
       def __hash__(self):
           MAX = int.MAX
           MASK = 2 * MAX + 1
           n = len(self)
           h = 1927868237 * (n + 1)
           h &= MASK
           for x in self:
               hx = hash(x)
               h ^= (hx ^ (hx << 16) ^ 89869747) * 3644798167
               h &= MASK
           h = h * 69069 + 907133923
           h &= MASK
           if h > MAX: h -= MASK + 1
           if h == -1: h = 590923713
           return h
   ```

2) **No quoted forward refs in annotations**  
   Use `Optional[TreeNode]`, **not** `Optional["TreeNode"]`. Same for `List[TreeNode]`, method return types, etc.

3) **Children must be lists**  
   Use `List[TreeNode]`. Don’t use tuples (Codon wants tuple sizes known at compile time; the tree is dynamic).

4) **Float everywhere**  
   Initialize numeric accumulators as `0.0`. Keep matrix type `List[List[float]]` to avoid Codon realizing functions on `int` first.

5) **Newick**  
   Implement a minimal parser and `Tree.from_newick()` as a **staticmethod**. Avoid a module-level `from_newick()` with the same name to prevent symbol shadowing.

6) **NJ final step**  
   Use a **3-child root** when `remaining == 3` to match Biotite’s reference, and sort the last three indices for determinism.

7) **Equality**  
   Don’t compare Newick strings. Compare all pairwise **topological** distances (exact) and **metric** distances (small epsilon).

8) **Codon NumPy I/O**  
   Avoid `np.loadtxt()` via Codon; use a pure file loader in `test_codon.py`.

---

## What the tests actually verify

- **`test_upgma`**: All pairwise leaf-to-leaf distances match a DendroUPGMA reference Newick within `1e-3`, and the **topology** (unweighted distances) matches exactly.
- **`test_neighbor_joining`**: NJ reconstruction equals a known reference tree (trifurcating root).
- **`test_distances`**: UPGMA is **ultrametric** (all leaves equidistant from root), and certain specific topological distances are correct.


