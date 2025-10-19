from typing import List, Optional

# For our tests we don't need a custom exception type
TreeError = RuntimeError
InvalidFileError = RuntimeError

# Codon has no frozenset and its set isn’t hashable by default.
# Provide a stable hash so we can compare unordered sets where needed.
@extend
class set:
    def __hash__(self):
        MAX = int.MAX
        MASK = 2 * MAX + 1
        n = len(self)
        h = 1927868237 * (n + 1)
        h &= MASK
        for x in self:
            hx = hash(x)
            h ^= (hx ^ (hx << 16) ^ 89869747) * 3644798167
            h &= MASK
        h = h * 69069 + 907133923
        h &= MASK
        if h > MAX:
            h -= MASK + 1
        if h == -1:
            h = 590923713
        return h


class TreeNode:
    _index: int                 # -1 for internal nodes, >=0 for leaves
    _distance: float            # branch length to parent
    _is_root: bool
    _parent: Optional[TreeNode]
    _children: List[TreeNode]

    def __init__(self,
                 children: Optional[List[TreeNode]] = None,
                 distances: Optional[List[float]] = None,
                 index: Optional[int] = None):
        self._is_root = False
        self._distance = 0.0
        self._parent = None
        if index is None:
            if children is None or distances is None:
                raise TypeError("Either index or (children, distances) must be set")
            if len(children) == 0:
                raise TreeError("Intermediate nodes need at least one child")
            if len(children) != len(distances):
                raise ValueError("children and distances length mismatch")
            # prevent the same object appearing twice
            for i in range(len(children)):
                for j in range(len(children)):
                    if i != j and children[i] is children[j]:
                        raise TreeError("Two child nodes cannot be the same object")
            self._index = -1
            self._children = [c for c in children]  # list copy; allow tuples/lists
            for k in range(len(self._children)):
                self._children[k]._set_parent(self, float(distances[k]))
        else:
            if index < 0:
                raise ValueError("Index cannot be negative")
            if children is not None or distances is not None:
                raise TypeError("Leaf node cannot have children/distances")
            self._index = int(index)
            self._children = []  # leaves have no children

    def _set_parent(self, parent: TreeNode, distance: float):
        if self._parent is not None or self._is_root:
            raise TreeError("Node already has a parent or is root")
        self._parent = parent
        self._distance = float(distance)

    @property
    def index(self) -> Optional[int]:
        return None if self._index == -1 else self._index

    @property
    def children(self) -> List[TreeNode]:
        return self._children

    @property
    def parent(self) -> Optional[TreeNode]:
        return self._parent

    @property
    def distance(self) -> float:
        # Biotite semantics: root reports 0.0
        return 0.0 if self._parent is None else self._distance

    def is_leaf(self) -> bool:
        return self._index != -1

    def is_root(self) -> bool:
        return self._is_root

    def as_root(self):
        if self._parent is not None:
            raise TreeError("Node has parent, cannot be a root node")
        self._is_root = True

    def get_leaves(self) -> List[TreeNode]:
        if self.is_leaf():
            return [self]
        out: List[TreeNode] = []
        for c in self._children:
            out.extend(c.get_leaves())
        return out

    def get_leaf_count(self) -> int:
        if self.is_leaf():
            return 1
        cnt = 0
        for c in self._children:
            cnt += c.get_leaf_count()
        return cnt

    def get_indices(self) -> List[int]:
        res: List[int] = []
        for lf in self.get_leaves():
            # lf.index cannot be None here
            res.append(lf.index)
        return res

    def _path_to_root(self) -> List[TreeNode]:
        path: List[TreeNode] = []
        cur: Optional[TreeNode] = self
        while cur is not None:
            path.append(cur)
            cur = cur._parent
        return path

    def lowest_common_ancestor(self, node: TreeNode) -> Optional[TreeNode]:
        a = self._path_to_root()
        b = node._path_to_root()
        i = len(a) - 1
        j = len(b) - 1
        lca: Optional[TreeNode] = None
        while i >= 0 and j >= 0:
            if a[i] is b[j]:
                lca = a[i]
                i -= 1
                j -= 1
            else:
                break
        return lca

    def distance_to(self, node: TreeNode, topological: bool = False) -> float:
        lca = self.lowest_common_ancestor(node)
        if lca is None:
            raise TreeError("Nodes do not share a common ancestor")
        dist = 0.0
        cur = self
        while cur is not lca:
            dist += 1.0 if topological else float(cur._distance)
            cur = cur._parent
        cur = node
        while cur is not lca:
            dist += 1.0 if topological else float(cur._distance)
            cur = cur._parent
        return dist

    def copy(self) -> TreeNode:
        if self.is_leaf():
            return TreeNode(index=self._index)
        clones = [c.copy() for c in self._children]
        dists = [float(c._distance) for c in self._children]
        return TreeNode(clones, dists)

    # --- Newick helpers -----------------------------------------------------

    def _label_for_leaf(self, labels: Optional[List[str]]):
        if labels is None:
            return str(self.index)
        idx = self.index
        if idx is None or idx < 0 or idx >= len(labels):
            raise ValueError("Leaf index out of range for provided labels")
        lab = labels[idx]
        # Keep labels Newick-safe (simple policy)
        if any(ch in lab for ch in [",", ":", ";", "(", ")", " "]):
            raise ValueError(f"Label contains illegal Newick characters: {lab!r}")
        return lab

    def to_newick(self,
                  labels: Optional[List[str]] = None,
                  include_distance: bool = True,
                  round_distance: Optional[int] = None) -> str:
        """
        Render this subtree to Newick.
        - labels: optional leaf-label list instead of numeric indices
        - include_distance: omit ':dist' when False
        - round_distance: round branch lengths to N decimals when set
        """
        def fmt_dist(x: float) -> str:
            if round_distance is not None:
                x = round(x, round_distance)
            if x == 0.0:
                x = 0.0
            return str(x)

        if self.is_leaf():
            label = self._label_for_leaf(labels)
            if include_distance:
                return f"{label}:{fmt_dist(self.distance)}"
            else:
                return label

        inner = ",".join(child.to_newick(labels, include_distance, round_distance)
                         for child in self._children)
        if include_distance:
            return f"({inner}):{fmt_dist(self.distance)}"
        else:
            return f"({inner})"

    # Pretty print (kept; not used by tests if to_newick is preferred)
    def __str__(self):
        if self.is_leaf():
            return f"{self._index}:{self._distance}"
        return "(" + ",".join(str(c) for c in self._children) + f"):{self._distance}"


class Tree:
    _root: TreeNode
    _leaves: List[TreeNode]

    def __init__(self, root: TreeNode):
        root.as_root()
        self._root = root
        leaves = root.get_leaves()
        n = len(leaves)
        if n == 0:
            raise TreeError("Tree must contain at least one leaf")
        # Pre-fill then place by index
        self._leaves = [leaves[0]] * n
        for lf in leaves:
            idx = lf.index
            if idx is None or idx < 0 or idx >= n:
                raise TreeError("The tree's indices are out of range")
            self._leaves[idx] = lf
        # ensure 0..n-1 mapping complete
        for i, lf in enumerate(self._leaves):
            if lf.index != i:
                raise TreeError("Leaf indices must be 0..n-1 without gaps")

    @property
    def root(self) -> TreeNode:
        return self._root

    @property
    def leaves(self) -> List[TreeNode]:
        return [lf for lf in self._leaves]

    @staticmethod
    def from_newick(s: str):
        # Parse: strip whitespace, drop trailing ';', parse to a node, wrap in Tree
        txt = "".join(c for c in s if not c.isspace())
        if len(txt) == 0:
            raise InvalidFileError("Newick string is empty")
        if txt.endswith(";"):
            txt = txt[:-1]
        node, _ = _parse_newick_node(txt)
        return Tree(node)

    def __len__(self) -> int:
        return len(self._leaves)

    def get_distance(self, i: int, j: int, topological: bool = False) -> float:
        return self._leaves[i].distance_to(self._leaves[j], topological)

    def copy(self):
        return Tree(self._root.copy())

    def to_newick(self,
                  labels: Optional[List[str]] = None,
                  include_distance: bool = True,
                  round_distance: Optional[int] = None) -> str:
        return self._root.to_newick(labels, include_distance, round_distance) + ";"

    # Robust, order-insensitive equality
    def __eq__(self, other) -> bool:
        if not isinstance(other, Tree):
            return False
        if len(self) != len(other):
            return False
        n = len(self)
        # Topology equality
        for i in range(n):
            for j in range(n):
                if self.get_distance(i, j, True) != other.get_distance(i, j, True):
                    return False
        # Metric equality with tiny tolerance
        EPS = 1e-9
        for i in range(n):
            for j in range(n):
                if abs(self.get_distance(i, j) - other.get_distance(i, j)) > EPS:
                    return False
        return True

    def __ne__(self, other) -> bool:
        return not self.__eq__(other)

    def __str__(self) -> str:
        return str(self._root) + ";"


# -------- Minimal Newick parser (sufficient for our tests) -------------------

def _parse_newick_node(s: str):
    if len(s) == 0:
        raise InvalidFileError("Empty node")
    if s[0] != "(":
        label = s
        dist = 0.0
        if ":" in s:
            parts = s.split(":", 1)
            label = parts[0]
            dist = float(parts[1])
        idx = int(label)
        return (TreeNode(index=idx), dist)

    level = 0
    close_idx = -1
    for i in range(len(s)):
        ch = s[i]
        if ch == "(":
            level += 1
        elif ch == ")":
            level -= 1
            if level == 0:
                close_idx = i
                break
    if close_idx == -1:
        raise InvalidFileError("Mismatched parentheses")

    inner = s[1:close_idx]
    tail = s[close_idx + 1:]
    dist = 0.0
    if len(tail) > 0:
        if tail[0] == ":":
            dist = float(tail[1:])
        elif ":" in tail:
            dist = float(tail.split(":", 1)[1])

    parts: List[str] = []
    buf = ""
    level = 0
    for ch in inner:
        if ch == "(":
            level += 1
            buf += ch
        elif ch == ")":
            level -= 1
            buf += ch
        elif ch == "," and level == 0:
            parts.append(buf)
            buf = ""
        else:
            buf += ch
    if len(buf) > 0:
        parts.append(buf)
    if len(parts) == 0:
        raise InvalidFileError("Intermediate node must have at least one child")

    children: List[TreeNode] = []
    dists: List[float] = []
    for p in parts:
        child, cd = _parse_newick_node(p)
        children.append(child)
        dists.append(cd)
    return (TreeNode(children, dists), dist)
