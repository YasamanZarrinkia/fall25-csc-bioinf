# AI Usage Log — Week 2

This document summarizes how AI assistance was used during the Week 2 TRViz (Codon) project. It groups prompts and intents by project phase and highlights recurring themes.

---

## Phase 1: Initial Setup and Understanding

**Prompts**
> "ok I have these sets of tests, I wanna know if this single test file covers all the test cases in the files I share later, by cases I mean if they all test exact same examples and set of datas"  
> "tell me what you changed and what you add and what you removed and why?"  
> "I want you to match my test file with expected decompositions, basically I want my test file to have exactly same data sets and examples and expected output"

**Intent / Outcome**
- Compared test coverage and ensured exact case/data parity across files.
- Explained diffs (added/removed/modified) and rationale.
- Synchronized expected outputs for perfect/imperfect/multi-motif cases and tie-breaks.

---

## Phase 2: Codon Type Errors and Fixes

**Prompts**
> "just give me the lines I need to change" (when getting Codon type errors)  
> "give me more details for test_decompose_dp_tie_break()"  
> "tell me how to create them on my repo" (asking about CI files)

**Intent / Outcome**
- Provided minimal diffs/line edits to satisfy Codon’s stricter typing.
- Clarified deterministic tie-break behavior and scoring combinations.
- Gave exact commands and file contents to create CI files in-repo.

---

## Phase 3: CI/CD Setup

**Prompts**
> "ok now since all my codes is in the repo my teacher give us an outline of how he marks it and we need to make a ci like below that runs our codes automatically" (shared CI template)  
> "are you running the test files in the evaluate? test_codon and test_python?"

**Intent / Outcome**
- Produced a working GitHub Actions workflow (robust Codon install + Python deps).
- Confirmed the evaluator runs both suites (`test_codon.py` and `test_python.py`).

---

## Phase 4: File Creation and Setup

**Prompts**
> "I already have week2 that it has the code folder in it with all the codon codes and the test files"  
> "ok write a report for me for this project, im sharing an example of my report for the week1 assignment" (shared Week 1 report)  
> "give me the .amd file" (meant ai.md)  
> "I don't want the ai.md for now, I want the report.md"  
> "give me all the prompts I used during this project"

**Intent / Outcome**
- Aligned paths to the existing `week2/code/` structure.
- Authored `REPORT.md` summarizing conversion scope and design choices.
- Generated `ai.md` (this file) and a consolidated list of user prompts.
- Ensured artifacts were downloadable for submission.

---

## Phase 5: Debugging and Error Resolution

**Implicit Prompts & Actions**
- Path resolution fixes (e.g., `PYTHONPATH`, colocating tests with code).
- Codon type error diagnostics and minimal-change patches.
- Test organization strategies (separate Codon/Python files with shared cases).
- Import structure guidance (avoid `trviz` in Codon; keep tests local to code).
- File structure organization for stable CI and reproducible local runs.

---

## Key Themes in Your Prompts

- **Precision:** Requests for exact matching outputs and line-specific changes.  
- **Debugging focus:** Targeted fixes for Codon and CI errors over high-level theory.  
- **Practical implementation:** Preference for concrete commands and complete file contents.  
- **CI/CD emphasis:** Strong focus on automation, reproducibility, and grading alignment.  
- **Teacher requirements:** Solutions mapped to marking criteria and assignment specs.

---

## Files Produced (Representative)

- `.github/workflows/week2.yml` — CI workflow running Codon + Python tests.  
- `week2/evaluate.sh` — Evaluator executing both suites with correct paths.  
- `week2/REPORT.md` — Technical report of the Week 2 work.  
- `week2/ai.md` — This AI usage summary.

---



