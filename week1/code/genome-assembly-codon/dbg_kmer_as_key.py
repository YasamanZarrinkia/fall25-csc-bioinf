from typing import Dict, Set, List, Optional

def reverse_complement(s: str) -> str:
    comp: Dict[str, str] = {'A':'T','T':'A','C':'G','G':'C'}
    out: List[str] = []
    for ch in reversed(s):
        out.append(comp.get(ch, ch))
    return ''.join(out)

class DBG_KmerKey:
    def __init__(self, k: int, data_list: List[List[str]]) -> None:
        self.k: int = k
        self.children: Dict[str, Set[str]] = {}
        self.parents: Dict[str, Set[str]] = {}
        self.counts: Dict[str, int] = {}
        for reads in data_list:
            for read in reads:
                self._add_read(read)

    def _add_read(self, read: str) -> None:
        k: int = self.k
        if len(read) < k:
            return
        prev: Optional[str] = None
        for i in range(len(read) - k + 1):
            kmer: str = read[i:i+k]
            self.counts[kmer] = self.counts.get(kmer, 0) + 1
            self.children.setdefault(kmer, set())
            self.parents.setdefault(kmer, set())
            if prev is not None and prev != kmer:
                if prev[1:] == kmer[:-1]:
                    self.children[prev].add(kmer)
                    self.parents[kmer].add(prev)
            prev = kmer

    def _walk_unitig_from(self, start: str) -> List[str]:
        path: List[str] = [start]
        cur: str = start
        while True:
            outs: List[str] = list(self.children.get(cur, set()))
            if len(outs) != 1:
                break
            nxt: str = outs[0]
            if len(self.parents.get(nxt, set())) != 1:
                break
            path.append(nxt)
            cur = nxt
        return path

    def _concat_path(self, path: List[str]) -> str:
        if not path:
            return ""
        s: str = path[0]
        for x in path[1:]:
            s += x[-1]
        return s

    def get_longest_contig(self) -> Optional[str]:
        best: str = ""
        starts: List[str] = [k for k in self.children.keys() if len(self.parents.get(k, set())) != 1 or len(self.children.get(k, set())) == 0]
        if not starts:
            starts = list(self.children.keys())
        for s in starts:
            path: List[str] = self._walk_unitig_from(s)
            seq: str = self._concat_path(path)
            if len(seq) > len(best):
                best = seq
        return best if best else None
