from typing import List, Tuple

def _join(path: str, name: str) -> str:
    # simple portable join without os.path
    if path.endswith("/") or path.endswith("\\"):
        return path + name
    return path + "/" + name

def read_fasta(path: str, name: str) -> List[str]:
    data: List[str] = []
    file_path: str = _join(path, name)
    with open(file_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line[0] != ">":
                data.append(line)
    return data

def read_data(path: str) -> Tuple[List[str], List[str], List[str]]:
    short1: List[str] = read_fasta(path, "short_1.fasta")
    short2: List[str] = read_fasta(path, "short_2.fasta")
    long1: List[str] = read_fasta(path, "long.fasta")
    return short1, short2, long1

