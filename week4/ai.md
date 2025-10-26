# AI Collaboration Log (Week 4)

**Model:** GPT-5 Thinking  
**Role requested:** “bioinformatic scientist who knows how to code in Python and Codon”

## Project summary
You asked me to implement four alignment algorithms (global, local, semi-global/fitting, and affine-gap global) in **separate Python and Codon implementations**, provide an `evaluate.sh` that prints a small timing table, run and explain results for `q2` vs `t2`, and set up an optional GitHub Actions CI workflow. We also debugged several environment and data issues together (multi-FASTA splitting, import paths, Codon typing, macOS `libomp.dylib`, etc.).

---

## Chronological prompt log (your requests)

1) **Initial assignment brief**  
   “You are a bioinformatic scientist… implement with Codon: global alignment; local alignment; semi-global (fitting) alignment; affine gap penalty global alignment. Also implement Python versions. Test with MT_human vs MT_orang and q1..q5 vs t1..t5 from ksw2. Use +3/−3/−2; affine −5 open, −1 extend. Provide `evaluate.sh` that prints:
   ```
   Method            Language    Runtime
   --------------------------------------
   global-mt_human   python      2000ms
   global-q1         python      2000ms
   local-q2          codon       1000ms
   ```”

2) **Separation of implementations**  
   “the python and codon implementation should be separate”

3) **Repo layout preference**  
   “create `week4` in my `fall25-csc-bioinf` repo and have all these under `week4/` (no extra outer folder).”

4) **Filesystem paths**  
   “the path to week4 is in Downloads while the fall-.. is in my Desktop”

5) **Download error (404) and fix**  
   You ran `./evaluate.sh` and got HTTP 404 during auto-download; asked how to fix (Option A: local copy, Option B: patch script).

6) **Clarification on Option A**  
   “should I change the evaluate.sh for optionA?”

7) **Filename mismatch**  
   You reported another 404 and later the underscore vs hyphen name issue (`MT_human.fa` vs `MT-human.fa`) and asked for a fix.

8) **Simplify evaluate & remove bench**  
   “this is my evaluate file, I also just wanna delete the bench since I already downloaded the datas” — asked for a minimal `evaluate.sh` that only times a few runs.

9) **Output formatting**  
   You showed mixed “score/time” and ms in the table; asked to make table show only `N ms`.

10) **Explain code & files**  
   “ok now I want you to tell me what you did with the codes, and walk me through each functions and each files”

11) **Data format realization (multi-FASTA)**  
   “q1.fa actually contains q1–q5; same for t1.fa (t1–t5)… there is no q2.fa etc.”

12) **Why Codon in one file vs Python split**  
   “so why are all the alignments for Codon in one file while Python is separated?”

13) **Force Codon-only for local-q2**  
   “for local-q2… I only want the Codon, no Python.”

14) **Installing/using Codon**  
   You hit Homebrew login prompt; asked for alternative install; then asked to avoid building a binary and use `codon run`.

15) **What’s ‘build a binary’ vs ‘codon run’?**  
   “what do you mean by building a binary or just use codon run”

16) **Border initialization in fitting**  
   “what do you mean by Fitting alignment’s borders: free first row (for t), penalized first column (for q)?”

17) **Show full results for q2 vs t2 (all methods, Py & Codon)**  
   “I want to just see the result … with both codon and python for all the methods.”

18) **Bug reports & fixes**  
   You shared Python SW TypeError → asked for fix; then Codon local typing error → asked for fix.

19) **Explain the outputs**  
   “ok now explain the output to me” — asked for interpretation of scores/alignments across methods.

20) **Why is Codon run so long**  
   “why my codon run is so long … `local-q2 codon 1392ms`”

21) **Share your current evaluate**  
   You pasted your `evaluate.sh` and asked how to reduce Codon timing (I suggested prefer binary if exists, else build, else run).

22) **CI (GitHub Actions) workflow**  
   “I need a .yml that runs my evaluate and is my CI; here’s a toy example.”

23) **Simplify CI**  
   “you didn’t need to do a lot; I already have the data.”

24) **Hard-require compiled binary**  
   “how can I do hard-require the compiled binary?”

25) **How to build the binary**  
   “how can I make the binary?”

26) **macOS `libomp.dylib` error**  
   You ran the binary and hit `Library not loaded: @loader_path/libomp.dylib` — asked for a self-contained fix.

27) **Self-contained evaluate**  
   “I prefer the self-contained approach, give me the updated evaluate file”

28) **Push week4 to repo**  
   “ok now I want to push my week4 folder to the repo”

29) **Create ai.md**  
   “actually give me an ai.md that has all the prompt I asked you to accomplish this project and it has your model”

---
