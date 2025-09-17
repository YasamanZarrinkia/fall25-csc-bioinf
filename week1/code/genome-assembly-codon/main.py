import sys
from typing import List, Optional
from dbg import DBG
from utils import read_data

def main() -> None:
    print("DEBUG: main start"); 
    argv: List[str] = sys.argv
    if len(argv) < 2:
        print("Usage: codon run main.py <data_dir>")
        sys.exit(1)

    data_dir: str = argv[1]
    short1, short2, long1 = read_data(data_dir)

    k: int = 25  # same k as Python
    dbg: DBG = DBG(k=k, data_list=[short1, short2, long1])

    # === DEBUG: quick graph stats ===
    node_count = len(dbg.nodes)
    edge_count = 0
    alive_nodes = 0
    for n in dbg.nodes:
        if n.alive:
            alive_nodes += 1
            edge_count += len(n.children)
    print("DEBUG graph:", "nodes=", node_count, "alive=", alive_nodes, "edges=", edge_count)
    # =================================

    out_path: str = data_dir + "/contig.fasta"
    with open(out_path, "w") as f:
        for i in range(20):
            c: Optional[str] = dbg.get_longest_contig()
            if c is None or len(c) == 0:
                break
            print(i, len(c))
            f.write(f">contig_{i}\n")
            f.write(f"{c}\n")

if __name__ == "__main__":
    main()

