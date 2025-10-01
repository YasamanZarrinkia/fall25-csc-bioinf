# Week 2 — TRViz (Codon) Port, Aligner/Encoder, and CI

## Overview
This week I ported the **DP decomposer** from TRViz to a Codon-friendly implementation and added lightweight, Codon-safe versions of the **motif encoder** and a **center-star–style aligner** (gap-padding only). I wired both Python and Codon test runners into CI. The Python tests rely on the `trviz` library; the Codon tests avoid it and exercise only my Codon-safe pipeline.

> **Note on tests:** I kept **two separate test files** (one for Codon, one for Python) although they cover the same cases, because Codon fails if a file imports `trviz`. Also, locating the tests inside the same directory as the code (`week2/code/`) keeps imports stable for Codon.

---

## Repository Layout (Week 2)
```
week2/
  ├─ code/
  │   ├─ decomposer.py        # DP decomposer + refine() (Codon-safe)
  │   ├─ main.py              # Worker wrapper, uses our encoder/aligner
  │   ├─ motif_aligner.py     # Center-star-like padding aligner (Codon-safe)
  │   ├─ motif_encoder.py     # Codon-safe encoder & score matrix helpers
  │   ├─ utils.py             # Counters, validation, edit distance, helpers
  │   ├─ test_codon.py        # Codon test runner (NO trviz import)
  │   └─ test_python.py       # Python test runner (uses trviz)
  └─ evaluate.sh              # CI runner: executes both test scripts
.github/workflows/week2.yml   # GitHub Actions (week2-ci)
```

---

## What Was Converted to Codon (1:1 behavior where feasible)

### `decomposer.py`
- `class Decomposer` (mode `"DP"` only).
- `decompose(sequence, motifs, …) -> List[str]`
  - Pure-Python DP (match/mismatch/insertion/deletion).
  - Deterministic tie-breaks to exactly match expected segmentations.
  - Strict validation (DNA alphabet; non-empty motif set).
  - No NumPy; fixed, Codon-friendly containers.
- Backtracking & reconstruction: equivalent behavior to TRViz DP in covered tests.
- `@staticmethod refine(decomposed_trs) -> List[List[str]]`
  - Majority/consensus per column; matches exact expected outputs.

### `motif_encoder.py`
- `get_motif_counter()`: frequency map (pure dict) for motif labels.
- **Labeling policy:**
  - Computes `private_motif_threshold` (`find_private_motif_threshold`) so total labels fit the alphabet.
  - Splits motifs into **normal** vs **private** (`_divide_motifs_into_normal_and_private`).
  - Assigns symbols from `INDEX_TO_CHR`, reserves `'?'` for private motifs.
- `encode(decomposed_vntrs, motif_map_file, label_count=None, auto=True)`
  - Writes a motif map (`write_motif_map`).
  - Returns concatenated **symbol strings** for each VNTR.
  - Persists `symbol_to_motif`, `motif_to_symbol`, and an empty `score_matrix` placeholder (Codon pipeline doesn’t use it, but API-compatible).
- `get_score_matrix(symbol_to_motif, …)`
  - Builds a simple substitution matrix using **Levenshtein distance**.
  - Gap penalties are included in the dict for downstream alignment heuristics.

### `motif_aligner.py`
- `class MotifAligner`
  - `align(sample_ids, encoded_vntrs, tool='center_star')`:
    - Implements center-star-like behavior using **gap padding** to the center’s length (first sequence as center).
    - Returns `(aligned_ids, aligned_seqs)`.
  - `_align_motifs_with_center_star`: simplified padding only; no DP multiple alignment.
  - `load_aligned_trs(aln_output)`: simple FASTA-like file parser fallback (optional).

### `utils.py`
- DNA checks: `is_valid_sequence(sequence)` (A/C/G/T only).
- Edit distance: `get_levenshtein_distance(s1, s2)` (pure Python).
- Indexing alphabet: `INDEX_TO_CHR`, `'?'` for private motifs.
- Padding and sorting:
  - `add_padding(encoded_trs)` to equalize lengths with `'-'`.
  - `sort(aligned_vntrs, sample_ids, symbol_to_motif, method='motif_count'|'name'|'manually')`.

---

## What Was Not Converted (and Why)

### TRViz HMM/CY modes
- Depend on `pomegranate`/NumPy; not Codon-compatible.
- A faithful port would be a separate project and may still diverge numerically.
- **Scope decision:** support **DP** mode only.

### Full TRViz visualizer / Matplotlib
- Plot generation stays on the Python side.
- Our CI tests **skip** alignment/plotting in Codon and perform basic padding alignment when needed.

### Importing `trviz` in Codon tests
- Codon fails on that import. Hence **separate test files**:
  - `test_codon.py` → only imports our Codon-safe modules.
  - `test_python.py` → imports TRViz and runs the same cases (plus encoder checks).

---

## Why Tests Aren’t in One File
I tried a combined test suite with conditional imports, but **Codon** still attempts to analyze the `trviz` import and fails.

Keeping two files guarantees:
- **Codon:** stable static analysis and imports (`from decomposer import Decomposer`, etc.).
- **Python:** full integration with TRViz where available.

Both files validate the **same sequences and expected outputs** for DP & refine, so behavior parity is maintained.

---

## Why Tests Live Next to the Code
Codon was “finicky” about relative imports and `sys.path` tricks.  
By placing tests in `week2/code/`, we can run both:

```bash
PYTHONPATH=week2/code codon run week2/code/test_codon.py
PYTHONPATH=week2/code python week2/code/test_python.py
```

with consistent imports in local runs and CI.

---

## CI & Tooling (week2-ci)

**Workflow:** `.github/workflows/week2.yml`
- Python 3.10 + `trviz` (for Python tests).
- Robust Codon installer that validates the tarball and adds it to PATH **in-step** and for later steps.
- Runs `week2/evaluate.sh` with `PYTHONPATH=week2/code`.

**Evaluator:** `week2/evaluate.sh`
- Checks for required files.
- Runs **Codon** tests first, then **Python** tests.
- Cleans small temporary artifacts (e.g., `*_motif_map.txt`).

---

## Behavioral Notes & Exactness
- **Decomposition parity:** For all perfect/imperfect cases and tie-break scenarios in the test suites, Codon and Python outputs are **exactly equal lists of motifs** (not just reconstruction equality).
- **Refine parity:** `refine()` produces the exact expected consensus sequences (tests assert full equality).
- **Alignment:** For Codon we intentionally use **padding** alignment to avoid DP multi-sequence alignment complexity and keep Codon dependencies pure. Python tests can still use TRViz’s encoder and then padding to maintain deterministic outputs.

---

## Known Limitations / Future Work
- **Center-star** here is a simplified alignment (gap-padding only). For real biological multiple alignment, integrate a DP-based progressive aligner in Python, or write a Codon-safe pairwise DP and progressive strategy (non-trivial).
- **HMM/CY** remain out of scope for Codon. If needed, re-implement with pure Python numerics and accept potential numeric differences.
- **Single-file tests:** feasible with heavy guards, but Codon’s static analysis still “sees” imports; keeping files split is more robust.

---

## Quick API Summary
- `Decomposer(mode="DP").decompose(sequence, motifs, **scores) -> List[str]`  
  Scores: `match_score`, `mismatch_score`, `insertion_score`, `deletion_score`.  
  Raises `ValueError` for invalid DNA or empty motifs.
- `Decomposer.refine(decomposed_trs) -> List[List[str]]`
- `MotifEncoder.encode(decomposed_vntrs, motif_map_file, label_count=None, auto=True) -> List[str]`  
  Side effects: writes motif map; sets `motif_to_symbol`, `symbol_to_motif`.
- `MotifAligner.align(sample_ids, encoded_vntrs, tool="center_star") -> (aligned_ids, aligned_seqs)`  
  (Padding alignment around the first sequence as center.)
- `utils`: `is_valid_sequence`, `get_levenshtein_distance`, `get_score_matrix`, `add_padding`, `sort`, constants.

---

## TL;DR
- **Converted to Codon:** DP decomposer, refine, motif encoder, simplified center-star aligner, and utilities — all pure Python and Codon-safe.  
- **Not converted:** HMM/CY modes, NumPy/pomegranate pieces, Matplotlib visualizer.  
- **Tests:** two files (Codon vs Python) for stability; they validate the same logic and cases.  
- **CI:** robust Codon install, Python deps, and a single evaluator running both paths.

PS: AI tools were used to help draft this report.