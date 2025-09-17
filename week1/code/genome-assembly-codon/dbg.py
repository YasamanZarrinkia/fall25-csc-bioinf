from typing import Dict, Set, List, Optional, Tuple

def reverse_complement(s: str) -> str:
    comp: Dict[str, str] = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    out: List[str] = []
    for ch in reversed(s):
        out.append(comp[ch] if ch in comp else ch)
    return ''.join(out)

class Node:
    kmer: str
    children: Set[int]
    count: int
    alive: bool
    visited: bool
    depth: int
    max_depth_child: Optional[int]
    
    def __init__(self, kmer: str):
        self.kmer = kmer
        self.children = set()
        self.count = 0
        self.alive = True
        self.visited = False
        self.depth = 0
        self.max_depth_child = None

class DBG:
    k: int
    nodes: List[Node]
    index: Dict[str, int]

    def __init__(self, k: int, data_list: List[List[str]]):
        self.k = k
        self.nodes = []
        self.index = {}
        self._build(data_list)

    def _build(self, data_list: List[List[str]]) -> None:
        for data in data_list:
            for original in data:
                rc = reverse_complement(original)
                # Exact match to original: range(len(original) - self.k - 1)
                for i in range(len(original) - self.k - 1):
                    self._add_arc(original[i:i+self.k], original[i+1:i+1+self.k])
                    self._add_arc(rc[i:i+self.k], rc[i+1:i+1+self.k])

    def _get_or_add(self, kmer: str) -> int:
        if kmer in self.index:
            return self.index[kmer]
        idx = len(self.nodes)
        self.index[kmer] = idx
        self.nodes.append(Node(kmer))
        return idx

    def _add_arc(self, kmer1: str, kmer2: str) -> None:
        idx1 = self._get_or_add(kmer1)
        idx2 = self._get_or_add(kmer2)
        self.nodes[idx1].count += 1
        self.nodes[idx2].count += 1
        if idx1 != idx2:
            self.nodes[idx1].children.add(idx2)

    def _reset(self) -> None:
        for node in self.nodes:
            if node.alive:
                node.visited = False
                node.depth = 0
                node.max_depth_child = None

    def _get_depth(self, idx: int) -> int:
        # Iterative DFS to replace recursion
        stack: List[Tuple[int, int]] = [(idx, 0)]  # (node_index, state)
        # state: 0 = not processed, 1 = children being processed
        
        while stack:
            node_idx, state = stack.pop()
            node = self.nodes[node_idx]
            
            if not node.alive:
                continue
                
            if state == 0:  # First visit
                node.visited = True
                # Push back with state=1 to process after children
                stack.append((node_idx, 1))
                
                # Get sorted children (like original)
                children = [c for c in node.children if self.nodes[c].alive]
                children.sort(key=lambda i: self.nodes[i].count, reverse=True)
                
                # Push children in reverse order to process in order
                for child_idx in reversed(children):
                    if not self.nodes[child_idx].visited:
                        stack.append((child_idx, 0))
            else:  # state == 1, process after children
                max_depth = 0
                max_child = None
                
                # Get sorted children again
                children = [c for c in node.children if self.nodes[c].alive]
                children.sort(key=lambda i: self.nodes[i].count, reverse=True)
                
                # Find child with maximum depth
                for child_idx in children:
                    child_node = self.nodes[child_idx]
                    if child_node.depth > max_depth:
                        max_depth = child_node.depth
                        max_child = child_idx
                
                # Set node's depth and max_depth_child
                node.depth = max_depth + 1
                node.max_depth_child = max_child
        
        return self.nodes[idx].depth

    def _get_longest_path(self) -> List[int]:
        max_depth = 0
        start_idx: Optional[int] = None
        
        for i, node in enumerate(self.nodes):
            if not node.alive:
                continue
            depth = self._get_depth(i)
            if depth > max_depth:
                max_depth = depth
                start_idx = i
        
        if start_idx is None:
            return []
        
        # Reconstruct path
        path: List[int] = []
        cur: Optional[int] = start_idx
        seen: Set[int] = set()
        
        while cur is not None and cur not in seen:
            if not self.nodes[cur].alive:
                break
            path.append(cur)
            seen.add(cur)
            cur = self.nodes[cur].max_depth_child
        
        return path

    def _concat_path(self, path: List[int]) -> str:
        if not path:
            return ""
        s = self.nodes[path[0]].kmer
        for nid in path[1:]:
            s += self.nodes[nid].kmer[-1]
        return s

    def _delete_path(self, path: List[int]) -> None:
        dead = set(path)
        for idx in path:
            node = self.nodes[idx]
            if node.kmer in self.index:
                del self.index[node.kmer]
            node.alive = False
        for n in self.nodes:
            if n.alive and n.children:
                n.children = n.children - dead

    def get_longest_contig(self) -> Optional[str]:
        self._reset()
        path = self._get_longest_path()
        if not path:
            return None
        contig = self._concat_path(path)
        self._delete_path(path)
        return contig