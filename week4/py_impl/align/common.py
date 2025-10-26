
from typing import Tuple, List

def read_fasta(path: str) -> str:
    seq = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line: continue
            if line.startswith(">"): 
                continue
            seq.append(line)
    return "".join(seq).upper()

def score(a: str, b: str, match: int, mismatch: int) -> int:
    return match if a == b else mismatch
