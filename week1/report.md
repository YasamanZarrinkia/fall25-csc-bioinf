# Week 1 Report

## Repository Layout
The repository is organized as follows:
```
.github/                  <- CI/CD configuration (GitHub Actions)
week1/
   code/
      genome-assembly/       <- Original Python source (provided)
      genome-assembly-codon/ <- Codon port of the source
   data/                     <- Input FASTA datasets (data1..data4)
   test/                     <- Output contigs written here by evaluate.sh
   ai.md                     <- AI tool identification and usage notes
   report.md                 <- This report
   evaluate.sh               <- Automated evaluation script
```

- `data/` contains unmodified input FASTA files (`short_1.fasta`, `short_2.fasta`, `long.fasta`)  
- `test/` stores evaluation outputs (`contig_py.fasta`, `contig_codon.fasta`), organized by dataset  
- `code/` contains both the original Python implementation and the Codon re-implementation  
- `evaluate.sh` runs both implementations, times them, and collects results  
- `ai.md` documents AI tools used during development and specific prompts/solutions  

---

## Setup and Environment
- **Python**: Runs locally or in Docker. Requires Python 3.11 and `matplotlib` (for plotting in the original code)  
- **Codon**: Can run locally (if installed) or in Docker. On Apple Silicon, Docker with `--platform=linux/amd64` is required  
- **Docker**: Ensures reproducibility across different machines and avoids local environment issues  

---

## Evaluation Script (`evaluate.sh`)
- Written in Bash (`#!/usr/bin/env bash`)  
- Begins with:
  ```bash
  set -euxo pipefail
  ```
  to ensure the script stops on errors, treats unset variables as failures, prints executed commands, and propagates failures in pipelines.

- The script:
  1. Locates the repository root (works from any directory)  
  2. Defines paths for `code/`, `data/`, `test/`, and `tmp/`  
  3. Supports running **locally** or inside **Docker**, controlled by environment variables:
     - `USE_DOCKER_PY=true` → force Python in Docker  
     - `USE_DOCKER_CODON=true` → force Codon in Docker  
     - `DOCKER_CPUS` / `DOCKER_MEM` → set Docker resources  
     - `DOCKER_PLATFORM_CODON=linux/amd64` → required on Apple Silicon  
  4. For each dataset (`data1..data4`), runs:
     - Python implementation → writes `contig_py.fasta`  
     - Codon implementation → writes `contig_codon.fasta`  
  5. Moves outputs into `week1/test/<dataset>/`  
  6. Computes N50 using a Python helper function  
  7. Prints a formatted summary table comparing runtime and N50 values

- Example run:
  ```bash
  bash week1/evaluate.sh
  ```

- Example with Docker and more resources:
  ```bash
  USE_DOCKER_CODON=true DOCKER_PLATFORM_CODON=linux/amd64 DOCKER_CPUS=6 DOCKER_MEM=16g bash week1/evaluate.sh
  ```

---

## Gotchas and Lessons Learned

- **Git workflow confusion:**  
  Initially attempted to clone the `genome-assembly` repository *inside* my `fall25-csc-bioinf` repo, creating complexity and confusion. With AI guidance, I adopted a cleaner approach: cloning locally, working on the code, and then pushing organized files into my main repository. This maintained a clean repository structure and simplified version control.

- **AI assistance:**  
  Heavily relied on AI tools throughout all assignment phases: repository setup, Python-to-Codon conversion, and evaluation automation. AI provided critical debugging assistance, especially for the Codon-specific recursion issue, and helped maintain organized project structure.

- **Dataset 4 challenge:**  
  Both implementations successfully processed `data4` when executed in Docker with sufficient resources. However, the initial Codon version consistently produced empty contig files for this dataset. Debugging revealed the root cause was a recursive depth-first search (DFS) implementation in `_get_depth()` that crashed with exit code 11 (stack overflow). This occurred because:
  - **Graph size:** De Bruijn graphs for genome assembly can be extremely large, leading to deep recursion
  - **Codon's stack limitations:** Codon has stricter stack constraints than standard Python, making it more susceptible to stack overflow with deep recursion

  The solution involved replacing the recursive DFS with an iterative stack-based approach that exactly replicated the original algorithm's behavior. After this fix:
  - Codon produced correct contigs with matching N50 values
  - Contig sequences were identical to Python's output
  - Only the ordering of contigs in the FASTA file sometimes differed (due to Codon's set/dict iteration order), which doesn't affect N50 calculations

- **Assembly Output and N50:**
  The assembler produces contigs in FASTA format (`contig.fasta`), with each entry containing an indexed header and the contig sequence. We maintained a fixed k-mer size of **k=25** for all runs, matching the original implementation.

  To evaluate assembly quality, we computed the **N50 statistic**:
  1. Sort all contigs by length (descending)  
  2. Calculate the total length of all contigs  
  3. Traverse the sorted list until reaching at least 50% of the total length  
  4. The length of the contig at this point is the **N50**  

  This was implemented as a helper function in `evaluate.sh` and applied to both Python and Codon outputs.

- **Performance observations:**
  As shown in the example output, Codon consistently outperformed Python in runtime, demonstrating the benefits of compiled execution for computationally intensive bioinformatics algorithms. The performance advantage was particularly noticeable for larger datasets like data4.

---

## Example Output

```
Dataset Language Runtime   N50
-------------------------------------------------------------------------------------------------------
data1   python   0:00:18   9990
data1   codon    0:00:08   9990
data2   python   0:00:33   9992
data2   codon    0:00:12   9992
data3   python   0:00:33   9824
data3   codon    0:00:13   9824
data4   python   0:08:58   159255
data4   codon    0:03:53   159255
```

---

## Key Takeaways

1. **Algorithm equivalence:** After resolving the recursion issue, the Codon implementation produces assembly results (N50) identical to the original Python version across all datasets.

2. **Performance advantage:** Codon demonstrates significant speed improvements (2-3x faster), especially valuable for larger datasets.

3. **Implementation considerations:** When porting algorithms to compiled environments like Codon, special attention must be paid to recursion depth and stack usage. Iterative implementations may be necessary for algorithms with deep recursion.

4. **Evaluation methodology:** N50 provides an effective metric for comparing assembly quality between implementations, as it focuses on contig lengths rather than ordering or implementation-specific details.

---

**PS:** AI tools were used to help draft this report.
