

## Prompts / questions I asked in chat (paraphrased)
GPT-5 Thinking

- “explain the following repo to me https://github.com/biotite-dev/biotite”  
- “what about https://www.biotite-python.org/latest/apidoc/biotite.sequence.phylo.html?”  
- “tell me what each of these three test files do, it’s testing related to the Phylo package” (and shared source)  
- “does it test all these you mentioned?”  
- “You are an expert Codon language programmer… (assignment details and constraints)… Modify the following functions so that they will be compatible with Codon rules… port biotite’s phylo package to Codon… run Python and Codon tests… print runtime table.”  
- “can you walk me through how to setup the repo with this layout? my repo name is `fall23-csc-bioinf` and instead of `phylo` I want `week3`.”  
- “Codon does not support unittest; change every test method to @test. Both codon test.py and python test.py should produce same results.” (toy example shared)  
- “I got: `ModuleNotFoundError: No module named 'week3'`”  
- “I don’t want a copy of codon codes in a .py file, the goal is to compare original Python codebase with converted Codon codes.”  
- “NumPy 2.0 error when importing biotite; what should I do?”  
- “Codon error: cannot import name 'join' from 'os.path' — are you sure it supports pathlib?”  
- “what should be the output when I run the test?”  
- “check if `REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))` … is codon friendly”  
- “still error with dirname — give me the full code”  
- “error: name 'globals' is not defined / `__file__` not defined / cannot import 'path' from 'sys'”  
- “you’re using `if HERE not in sys.path` again”  
- “unknown arg `-I` in codon run; how to set include paths?”  
- “unexpected dedent / expected type expression / syntax errors in `tree.codon`”  
- “im tired of getting all these errors, go find how properly write codon friendly code”  
- “toy example for detecting Codon vs Python (memcpy on str); use that style”  
- “errors with exceptions must derive from BaseException / Static[Exception]”  
- “same error; give me the full code”  
- “still errors in `__init__` of TreeNode; children/distance mismatch; unexpected dedent”  
- “hashing/equality issues with sets and TreeNode; errors with '!='”  
- “CError: dlopen(libpython.dylib) — remove Python bridge”  
- “test_codon.py type mismatch: List[List[float]] vs List[List[int]]; fix”  
- “Runtime table shows python and codon times; codon initially 0ms — then fixed”  
- “why don’t you use the .txt files in the test folder?”  
- “here is the original test file; please convert our tests to those 3 items”  
- “TypeError in pytest_approx mock; fix”  
- “Hints: load distances.txt with python bridge or pure loader; use float64; set hashing for set; use RandomError alias… (I shared instructor hints)”  
- “this is my `tree.codon`; this is `upgma.codon`; this is `neighbor_joining.codon`”  
- “give me the full code for all these updated files”  
- “cannot import typing.Static — fix”  
- “share a downloadable report.md (steps & gotchas)”  
- “what I meant was: share the commands/prompts I asked you in chat so I can give them to my instructor.”

---



