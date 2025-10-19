# Only used when week3 is imported under Codon; Python tests import Biotite.
try:
    __codon__  # type: ignore
    from .tree import Tree, TreeNode, TreeError, InvalidFileError
    from .upgma import upgma
    from .neighbor_joining import neighbor_joining
except NameError:
    pass
