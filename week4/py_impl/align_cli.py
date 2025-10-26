import argparse, time
from .align.common import read_fasta
from .align.nw import global_alignment
from .align.sw import local_alignment
from .align.fit import semi_global_fitting
from .align.affine import affine_global

def main():
    p = argparse.ArgumentParser(description="Python DP aligners")
    p.add_argument("--method", choices=["global","local","fitting","affine"], required=True)
    p.add_argument("--query", required=True)
    p.add_argument("--target", required=True)
    p.add_argument("--match", type=int, default=3)
    p.add_argument("--mismatch", type=int, default=-3)
    p.add_argument("--gap", type=int, default=-2)
    p.add_argument("--gap_open", type=int, default=-5)
    p.add_argument("--gap_ext", type=int, default=-1)
    args = p.parse_args()
    q = read_fasta(args.query)
    t = read_fasta(args.target)
    t0 = time.time()
    if args.method=="global":
        sc, qa, ta = global_alignment(q,t,args.match,args.mismatch,args.gap)
    elif args.method=="local":
        sc, qa, ta = local_alignment(q,t,args.match,args.mismatch,args.gap)
    elif args.method=="fitting":
        sc, qa, ta = semi_global_fitting(q,t,args.match,args.mismatch,args.gap)
    else:
        sc, qa, ta = affine_global(q,t,args.match,args.mismatch,args.gap_open,args.gap_ext)
    ms = int((time.time()-t0)*1000)
    print(f"score={sc} time_ms={ms}")
    if len(qa) <= 200:
        print(qa); print(ta)
    else:
        print("(alignment suppressed)")
if __name__ == "__main__":
    main()
