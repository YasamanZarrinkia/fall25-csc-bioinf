
from typing import List, Dict, Tuple, Iterable
import string
from collections import Counter

# --------- Constants / Symbols ----------
LOWERCASE_LETTERS = string.ascii_lowercase
UPPERCASE_LETTERS = string.ascii_uppercase
DIGITS = string.digits

skipping_characters = ['(', '=', '<', '>', '?', '-']
PRIVATE_MOTIF_LABEL = '?'
INDEX_TO_CHR = list(LOWERCASE_LETTERS) + list(UPPERCASE_LETTERS) + list(DIGITS)
INDEX_TO_CHR.extend([chr(x) for x in range(33, 127) if chr(x) not in skipping_characters and chr(x) not in INDEX_TO_CHR])

DNA_CHARACTERS = {'A','C','G','T'}

# --------- FASTA I/O (no Biopython) ----------
def get_sample_and_sequence_from_fasta(fasta_file: str) -> Tuple[List[str], List[str]]:
    headers: List[str] = []
    seqs: List[str] = []
    curr = []
    curr_h = None
    with open(fasta_file, "r") as f:
        for line in f:
            line=line.strip()
            if not line: continue
            if line.startswith('>'):
                if curr_h is not None:
                    seqs.append(''.join(curr).upper())
                    curr = []
                curr_h = line[1:].strip()
                headers.append(curr_h)
            else:
                curr.append(line)
        if curr_h is not None:
            seqs.append(''.join(curr).upper())
    return headers, seqs

# --------- Basic validations ----------
def is_valid_sequence(seq: str) -> bool:
    s = seq.upper()
    return all(c in DNA_CHARACTERS for c in s)

# --------- Motif helpers ----------
def get_motif_counter(decomposed_vntrs: List[List[str]]) -> Counter:
    c = Counter()
    for d in decomposed_vntrs:
        for m in d:
            c[m] += 1
    return c

# --------- Distances ----------
def levenshtein(a: str, b: str) -> int:
    n, m = len(a), len(b)
    dp = [list(range(m+1))] + [[i]+[0]*m for i in range(1,n+1)]
    for i in range(1,n+1):
        ai=a[i-1]
        row=dp[i]
        prow=dp[i-1]
        for j in range(1,m+1):
            cost = 0 if ai==b[j-1] else 1
            row[j] = min(prow[j]+1, row[j-1]+1, prow[j-1]+cost)
    return dp[n][m]

# --------- Score matrix over symbols ----------
def get_score_matrix(symbol_to_motif: Dict[str,str],
                     match_score: float = 2.0,
                     mild_mismatch: float = -1.0,
                     harsh_mismatch: float = -2.0,
                     gap_open: float = 1.5,
                     gap_extension: float = 0.6) -> Dict:
    symbols = list(symbol_to_motif.keys())
    score = {s1:{} for s1 in symbols + ['?'] if s1 not in {}}
    # ensure '?' present
    if '?' not in score:
        score['?'] = {}
    for s1 in score:
        for s2 in score:
            if s1 == s2:
                score[s1][s2] = match_score
            else:
                if s1=='?' or s2=='?':
                    score[s1][s2] = harsh_mismatch
                else:
                    m1 = symbol_to_motif[s1]
                    m2 = symbol_to_motif[s2]
                    d = levenshtein(m1, m2)
                    cutoff = 1 + max(len(m1), len(m2))//30
                    score[s1][s2] = mild_mismatch if d <= cutoff else harsh_mismatch
    score['gap_open'] = gap_open
    score['gap_extension'] = gap_extension
    return score

# --------- Sorting (minimal) ----------
def sort(sample_ids: List[str], aligned_vntrs: List[str], method: str = 'name', sample_order_file: str = None):
    idxs = list(range(len(sample_ids)))
    if method == 'name':
        idxs.sort(key=lambda i: str(sample_ids[i]))
    elif method == 'motif_count':
        idxs.sort(key=lambda i: (-sum(1 for c in aligned_vntrs[i] if c!='-'), str(sample_ids[i])))
    else:
        # no-op for unsupported methods in Codon build
        pass
    return [sample_ids[i] for i in idxs], [aligned_vntrs[i] for i in idxs]

# --------- Minor layout helpers (stubs) ----------
def add_padding(aligned_vntrs: List[str], pad_left: int = 0, pad_right: int = 0) -> List[str]:
    left = '-' * pad_left
    right = '-' * pad_right
    return [left + s + right for s in aligned_vntrs]

def get_motif_marks(aligned_vntrs: List[str]) -> List[int]:
    # returns column indices with non-gap content (for simple ticks)
    if not aligned_vntrs: return []
    L = len(aligned_vntrs[0])
    marks = [i for i in range(L) if any(r[i] != '-' for r in aligned_vntrs)]
    return marks

