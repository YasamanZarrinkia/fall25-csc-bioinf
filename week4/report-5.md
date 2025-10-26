# Week 4 — Edit Distance & Dynamic Programming Alignments (Python & Codon)

This report documents **what was built**, **how to run it**, **why each algorithm behaves the way it does**, and **every gotcha we hit**.

---

## 1) Goal & Requirements

Implement the following alignment algorithms in **both Python and Codon**, then evaluate on the **ksw2** test sequences using the scoring:
- Linear-gap methods: **+3** match, **−3** mismatch, **−2** gap
- Affine-gap (Gotoh): **+3** match, **−3** mismatch, **−5** gap-open, **−1** gap-extend

Algorithms:
1. Global alignment (Needleman–Wunsch; linear gap)
2. Local alignment (Smith–Waterman; linear gap)
3. Semi-global / Fitting alignment (linear gap)
4. Global alignment with affine gaps (Gotoh; 3 matrices)


---

## 2) Final Project Layout

```
week4/
  py_impl/
    __init__.py
    align/
      common.py
      nw.py
      sw.py
      fit.py
      affine.py
    align_cli.py
  codon_impl/
    align_codon.py      # Codon-only, standalone CLI with all 4 methods
  data/                 # FASTA inputs 
  evaluate.sh           # prints the 3-row timing table
```

**Why Python is split but Codon is single-file:** For Codon I favored a single file to keep the build/run dead-simple and avoid import path pitfalls; Python is idiomatic and easier to read in modules. 
---

## 3) Data Inputs (ksw2)

The ksw2 test inputs are **multi-FASTA** files (`q1.fa` and `t1.fa`) that contain **q1…q5** and **t1…t5** entries. For convenience, I **split** them into per-entry files so the CLIs can be called as:

```
data/q2.fa  vs  data/t2.fa
```

---

## 4) Algorithms — How They’re Implemented

### 4.1 Global Alignment (Needleman–Wunsch, linear gap)
- DP table `dp[i][j]`: best score aligning `a[:i]` to `b[:j]`.
- Initialization (penalize borders):
  - `dp[i][0] = i * gap` (gaps in `b`)
  - `dp[0][j] = j * gap` (gaps in `a`)
- Recurrence at cell `(i,j)`:
  - `diag = dp[i-1][j-1] + score(a[i-1], b[j-1])`
  - `up   = dp[i-1][j] + gap`     (gap in `b`)
  - `left = dp[i][j-1] + gap`     (gap in `a`)
  - `dp[i][j] = max(diag, up, left)`
- Traceback from `(n,m)`. Complexity: `O(nm)`.

### 4.2 Local Alignment (Smith–Waterman, linear gap)
- Allows the alignment to start/end anywhere inside the sequences.
- Same recurrence but includes **0 reset**: `dp[i][j] = max(0, diag, up, left)`.
- Track the **best cell** `(bi, bj)`; traceback stops at the first zero. Complexity `O(nm)`.

### 4.3 Semi-global / Fitting Alignment (linear gap)
- Fit **entire `q`** inside a **substring of `t`**.
- Border decisions (“fitting borders”):
  - **Top row (i = 0) free:** `dp[0][j] = 0` (you can start in `t` anywhere)
  - **First column penalized:** `dp[i][0] = i*gap` (must consume all of `q`)
- End: take the best over the **last row** `max_j dp[len(q)][j]`.
- Traceback from `(len(q), j_best)`; any remaining `t` overhang is **free**.

### 4.4 Global Alignment with Affine Gaps (Gotoh)
- Three matrices:
  - `M[i][j]`: ends with a **match/mismatch**
  - `X[i][j]`: ends with a **gap in b** (up move)
  - `Y[i][j]`: ends with a **gap in a** (left move)
- Transitions:
  - `X[i][j] = max(M[i-1][j] + gap_open,  X[i-1][j] + gap_ext)`
  - `Y[i][j] = max(M[i][j-1] + gap_open,  Y[i][j-1] + gap_ext)`
  - `M[i][j] = max(M[i-1][j-1], X[i-1][j-1], Y[i-1][j-1]) + score`
- Initialization:
  - `M[0][0] = 0`, first row/col set to reflect gap-open then extend.
- End: `max(M[n][m], X[n][m], Y[n][m])`. Traceback across matrices.

---

## 5) CLIs

### 5.1 Python CLI
```
python -m py_impl.align_cli --method {global|local|fitting|affine}   --query data/MT-human.fa --target data/MT-orang.fa   [--match 3 --mismatch -3 --gap -2 --gap_open -5 --gap_ext -1]
```
It prints `score=... time_ms=...` and the alignment (suppressed only for very long outputs).

### 5.2 Codon CLI
- Build & run (binary):  
  `codon build -release -o bin/align_codon codon_impl/align_codon.py && ./bin/align_codon ...`
- Or run-in-place (no binary):  
  `codon run -release codon_impl/align_codon.py --method local --query data/q2.fa --target data/t2.fa`

> **Binary vs `codon run`:** Building creates a fast standalone executable; `codon run` compiles-then-executes on the fly (no artifact), great for iteration.

---

## 6) The small timing table (`evaluate.sh`)

The requested 3-row output:
```
Method            Language    Runtime
--------------------------------------
global-mt_human   python      <N>ms
global-q1         python      <N>ms
local-q2          codon       <N>ms
```

---

## 8) Gotchas I Hit (and Fixes)

1. **Network/RAW 404s:**  
   `raw.githubusercontent.com` can intermittently 404. I eliminated downloads entirely by copying FASTAs locally and, alternatively, provided fallback URLs.

2. **Multi-FASTA inputs:**  
   `q1.fa`/`t1.fa` contain multiple entries (`q1…q5`, `t1…t5`). We added a small **splitter** to generate `q2.fa`, `t2.fa`, etc., to keep CLIs simple.

3. **Python imports:**  
   The CLI originally imported from `py_impl.common`; fixed by making `py_impl` a package and importing from `.align.*`.

---

### TL;DR
 Now have clean, separate Python and Codon implementations for four DP alignment methods, correct scoring, a small reproducible timing table. All the common pitfalls (filenames, multi-FASTA, imports, strict typing, output formatting) have been handled and documented here.
