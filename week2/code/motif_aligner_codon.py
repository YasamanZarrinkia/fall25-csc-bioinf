
import subprocess
import os
from typing import Tuple, List, Dict

class MotifAligner:
    def align(self,
              sample_ids: List[str],
              encoded_vntrs: List[str],
              vid: str = None,
              score_matrix: Dict = None,
              output_dir: str = "./",
              tool: str = "mafft",
              ) -> Tuple[List, List]:
        motif_aligner = self._get_motif_aligner(tool)
        return motif_aligner(sample_ids, encoded_vntrs, vid, score_matrix, output_dir)

    def _get_motif_aligner(self, tool):
        if tool == 'mafft':
            return self._align_motifs_with_mafft
        elif tool == 'star':
            return self._align_motifs_with_star
        else:
            # default to mafft, then fallback to star
            return self._align_motifs_with_mafft

    @staticmethod
    def _write_fasta(sample_ids: List[str], labeled_vntrs: List[str]) -> str:
        from tempfile import mkstemp
        fd, path = mkstemp(prefix=\"trviz_\", suffix=\".fa\")
        os.close(fd)
        with open(path, \"w\") as f:
            for sid, lab in zip(sample_ids, labeled_vntrs):
                f.write(f\">{sid}\\n{lab}\\n\")
        return path

    @staticmethod
    def _write_mafft_matrix(score_matrix: Dict) -> str:
        # Writes a MAFFT-compatible matrix for our alphabet (minus gaps). Gaps are set via flags.
        from tempfile import mkstemp
        fd, path = mkstemp(prefix=\"trviz_\", suffix=\".mat\")
        os.close(fd)
        symbols = [s for s in score_matrix.keys() if len(s)==1 and s.isprintable() and s not in ('gap_open','gap_extension')]
        with open(path, \"w\") as f:
            f.write(\"  \" + \" \".join(symbols) + \"\\n\")
            for s1 in symbols:
                row = [s1] + [str(int(score_matrix[s1][s2])) for s2 in symbols]
                f.write(\" \".join(row) + \"\\n\")
        return path, symbols

    def _align_motifs_with_mafft(self, sample_ids, labeled_vntrs, vid, score_matrix, output_dir):
        # Try to call mafft; if not present, fallback to star aligner
        try:
            in_fa = self._write_fasta(sample_ids, labeled_vntrs)
            mat_path, symbols = self._write_mafft_matrix(score_matrix if score_matrix else {'?': {'?':2}})
            cmd = [\"mafft\", \"--text\", \"--matrix\", mat_path, in_fa]
            res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False, text=True)
            if res.returncode != 0:
                # fallback
                return self._align_motifs_with_star(sample_ids, labeled_vntrs, vid, score_matrix, output_dir)
            # parse alignment
            aligned_ids = []
            aligned = []
            sid = None
            seq = []
            for line in res.stdout.splitlines():
                line=line.strip()
                if not line: continue
                if line.startswith('>'):
                    if sid is not None:
                        aligned_ids.append(sid); aligned.append(''.join(seq))
                        seq=[]
                    sid = line[1:].strip()
                else:
                    seq.append(line)
            if sid is not None:
                aligned_ids.append(sid); aligned.append(''.join(seq))
            # preserve the original order
            id_to_row = {sid: row for sid, row in zip(aligned_ids, aligned)}
            aligned_out = [id_to_row[s] for s in sample_ids]
            return sample_ids, aligned_out
        except Exception:
            return self._align_motifs_with_star(sample_ids, labeled_vntrs, vid, score_matrix, output_dir)

    # --- simple star alignment fallback (Needleman-Wunsch) ---
    @staticmethod
    def _nw(a: str, b: str, score: Dict, gap_open: float = 1.5) -> tuple:
        gap = -gap_open
        n, m = len(a), len(b)
        dp = [[0.0]*(m+1) for _ in range(n+1)]
        bt = [[(0,0)]*(m+1) for _ in range(n+1)]
        for i in range(1, n+1):
            dp[i][0] = dp[i-1][0] + gap; bt[i][0] = (i-1,0)
        for j in range(1, m+1):
            dp[0][j] = dp[0][j-1] + gap; bt[0][j] = (0,j-1)
        for i in range(1, n+1):
            for j in range(1, m+1):
                sdiag = dp[i-1][j-1] + score[a[i-1]][b[j-1]]
                sup   = dp[i-1][j] + gap
                sleft = dp[i][j-1] + gap
                best = sdiag; prev=(i-1,j-1)
                if sup > best: best, prev = sup, (i-1,j)
                if sleft > best: best, prev = sleft, (i,j-1)
                dp[i][j] = best; bt[i][j] = prev
        i,j = n,m
        A,B=[],[]
        while i>0 or j>0:
            pi,pj = bt[i][j]
            if pi==i-1 and pj==j-1:
                A.append(a[i-1]); B.append(b[j-1])
            elif pi==i-1 and pj==j:
                A.append(a[i-1]); B.append('-')
            else:
                A.append('-'); B.append(b[j-1])
            i,j = pi,pj
        return ''.join(reversed(A)), ''.join(reversed(B)), dp[n][m]

    def _align_motifs_with_star(self, sample_ids, labeled_vntrs, vid, score_matrix, output_dir):
        # pick naive center: the first string
        center = labeled_vntrs[0]
        aligned = [center]
        ids = [sample_ids[0]]
        for sid, s in zip(sample_ids[1:], labeled_vntrs[1:]):
            A,B,_ = self._nw(aligned[0], s, score_matrix if score_matrix else {'?': {'?':2}}, gap_open=score_matrix.get('gap_open',1.5) if score_matrix else 1.5)
            # pad previous aligned to A's gaps
            def expand(template, s2):
                it = iter(s2)
                out=[]
                for ch in template:
                    if ch=='-':
                        out.append('-')
                    else:
                        out.append(next(it))
                out.extend(list(it))
                return ''.join(out)
            aligned = [expand(A, aligned[0].replace('-',''))]
            aligned.append(B)
            ids = [ids[0], sid]
        return ids, aligned

