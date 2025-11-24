# AI Usage Statement – Week 6 (Single-cell RNA-seq)

For this assignment I used **ChatGPT (OpenAI, GPT-5.1 Thinking via ChatGPT Plus)** as a coding/debugging assistant and for help with explanation and wording. All code, commands, and analysis steps were **run, checked, and adapted by me** on my own machine and in my GitHub CI. The AI did not have access to my data or to my GitHub repository; I copied code and error messages manually into the chat and then copied the relevant suggestions back into my notebook and workflow files.

Below I summarize **how** I used AI and give **representative questions** I asked while completing this project.

---

## 1. Environment, Jupyter, and basic tooling

I used AI to help with environment setup and basic Jupyter usage, especially when I got stuck on small things that are easy to forget.

The assistant answered with terminal commands (e.g., `jupyter notebook`, `cd week6`, `mkdir`, `git rm -r` etc.), and short explanations. I executed these myself and adjusted paths as needed.

---

## 2. Single-cell pipeline design & debugging (Scanpy, CellTypist, QC)

I used AI while building and debugging the single-cell RNA-seq pipeline in the notebook (Scanpy, CellTypist, QC plots, etc.). I already followed the **Single-cell Best Practices Book**, but the assistant helped me interpret errors and adjust parameters.

Representative questions:

- When filtering removed all cells:

  > “Filtering cells and genes…  
  > Before filtering: 139 cells, 60 genes  
  > filtered out 139 cells that have less than 200 genes expressed  
  > filtered out 60 genes that are detected in less than 3 cells  
  > After filtering: 0 cells, 0 genes  
  > That’s not good, now the visualization is not complete — I’m confused. What should I do?”

- When Scanpy QC crashed:

  > “Calculating QC metrics… `IndexError: Positions outside range of features.`  
  > What does this mean and how can I fix it?”

- For highly variable genes:

  > “Finding highly variable genes… `ImportError: Please install skmisc package` — do I really need this for a small toy dataset, or can I change the flavor/parameters?”

- For PCA:

  > “PCA error: `ValueError: n_components=50 must be between 1 and min(n_samples, n_features)=30` — what should I set for `n_comps` in this toy example?”

- For CellTypist:

  > “CellTypist error: `ValueError: No features overlap with the model. Please provide gene symbols`.  
  > I have 60 genes and Ensembl-like IDs. How can I map them to gene symbols so CellTypist works?”

The assistant suggested:

- Relaxing filtering thresholds for the toy dataset (e.g., requiring fewer genes per cell).
- Using a different `flavor` for `sc.pp.highly_variable_genes` that does not require `skmisc` or adjusting parameters.
- Reducing `n_comps` in PCA to a safe value (e.g., 20–30) based on `min(n_cells, n_genes)`.
- Writing helper code to map Ensembl IDs to gene symbols using my GTF (and then using those gene symbols in CellTypist).

I implemented these suggestions, reran the notebook, inspected UMAPs, QC plots, and verified that CellTypist successfully produced labels (or at least ran without errors on the toy data).

---

## 3. Alevin-fry, whitelist, and Salmon / alevin issues

I used AI quite heavily to debug the **alevin-fry** and **Salmon alevin** steps, since the error messages are quite low-level.

Representative questions:

- When using the external whitelist:

  > “`alevin-fry generate-permit-list` crashes with  
  > `assertion 'left == right' failed: found barcodes of different lengths 42 and 75`.  
  > My local `3M-february-2018.txt` looks fine (all lines length 16). What is going on?”

- When I tried knee-distance:

  > “I changed the command to use `--knee-distance` (no whitelist), but in CI it says  
  > `Error: could not open input rad file` even though it works on my machine. Why is this happening?”

- For Salmon index building:

  > “Salmon index fails in CI with `terminate called without an active exception` and mentions `Filter size` and `k-mers`. What might cause this on the toy splici FASTA, and can I instead pre-compute the index locally and commit it?”

- When the notebook couldn’t find `quants_mat.mtx`:

  > “Loading alevin-fry output… FileNotFoundError: `data/alevin_fry_quant/alevin/quants_mat.mtx` not found in CI. It exists on my laptop. How do I make CI see it?”

The assistant helped me:

- Understand that the 42 vs 75 barcode length error likely came from trying to use the big 3M whitelist in a way that alevin-fry didn’t like, and that using the **knee-distance** mode or a properly formatted whitelist would be safer.
- Realize that CI had **no RAD file** because Salmon hadn’t run (or the index failed), so the RAD-dependent steps crashed.
- Decide on a strategy that is consistent with the assignment + teacher’s note:
  - **Compute the Salmon index and alevin + alevin-fry outputs locally** with the toy data.
  - Check file sizes to keep them below GitHub’s 100 MB limit.
  - Commit **only** those small index + quant outputs (not the huge reference genome FASTA) so that CI can **skip** heavy alignment and just load pre-computed outputs.
- Write idempotent bash cells like:
  - “If `data/salmon_index/versionInfo.json` exists, skip building the index.”
  - “If `data/alevin_output/map.rad` exists, skip running Salmon.”
  - “If `data/alevin_fry_quant/alevin/quants_mat.mtx` exists, skip running alevin-fry.”

This kept the notebook reproducible for CI while avoiding the fragile steps on GitHub runners.

---

## 4. Git, GitHub CI, YAML, and workflow design

I also used AI to help design and debug the **GitHub Actions CI** (`.github/workflows/week6.yml`), following the instructor’s requirement that `jupyter execute week6.ipynb` (or equivalent) should reproduce everything.

Representative questions:

- “This is my current CI YAML (pasted). What do you think? I’m not sure if I should use micromamba or conda-incubator.”
- “I get `Error: micromamba-version must be either latest or a version matching 1.2.3-0`. What value should I use for `micromamba-version`?”
- “I have this toy example CI with Codon. Can you give me a proper CI that sets up an environment with salmon, alevin-fry, scanpy, etc., and then runs `week6/week6.ipynb`?”
- “I get `ModuleNotFoundError: No module named 'scanpy'` in CI even though I installed things. What did I do wrong in the YAML?”
- “GitHub says: `Invalid workflow file .github/workflows/week6.yml#L47 – You have an error in your yaml syntax on line 47`. Can you help me fix the indentation?”
- “My commit failed because `genome_fixed.fa` is 175 MB and exceeds GitHub’s 100 MB limit. How can I undo that commit and avoid pushing huge files?”

The assistant suggested:

- A micromamba-based CI workflow that:
  - Uses `mamba-org/setup-micromamba@v2`.
  - Installs `salmon`, `alevin-fry`, `scanpy`, `anndata`, `celltypist`, `pyroe`, `jupyter`, and basic scientific Python stack.
  - Runs:

    ```bash
    cd week6
    jupyter nbconvert --to notebook --execute week6.ipynb --output week6_executed.ipynb
    ```

  - Uploads the executed notebook as an artifact.

- Correct YAML indentation and shell syntax.
- Using **guards** so that the notebook can safely re-run on CI using pre-computed data.
- Avoiding very large files by:
  - Committing only the Salmon index and quant outputs that are small enough.
  - Removing accidentally committed large files from history when needed (using `git reset` / `git rm` etc., which I ran myself).

---

## 5. Report text and markdown explanation

Finally, I used AI to help me **word** parts of the notebook markdown and write a concise description of:

- The overall pipeline (from FASTQ → alevin → alevin-fry → AnnData → clustering → annotation).
- Interpretation of QC plots, PCA, UMAP, and marker gene plots.
- A short write-up on the CellTypist results and how they relate to the clusters.

Representative questions:

- “Just give me a complete report so I can add it in the markdown: explain what I did with QC, filtering, HVG selection, PCA, UMAP, clustering and CellTypist on this toy dataset.”
- “Can you summarize the changes we made (knee-distance, precomputed outputs, etc.) so I can explain in the notebook why CI still reproduces the analysis?”

I edited this text to match my own style and to accurately reflect what I actually ran.

---

## 6. My responsibility

Even though I relied on ChatGPT for help with:

- Environment setup,
- Error interpretation,
- Example YAML configs,
- and wording in markdown,

I:

- Ran all commands myself and verified that they work.
- Manually checked the CI workflows and logs on GitHub.
- Reviewed and edited all code that ended up in `week6/week6.ipynb` and `.github/workflows/week6.yml`.
- Ensured that the final notebook **really does** implement:
  1. Alignment and quantification (via Salmon alevin + alevin-fry, with pre-computed outputs committed),
  2. Clustering (Leiden),
  3. Automatic annotation (CellTypist),
  4. Visualization and reporting.

