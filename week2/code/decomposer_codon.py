from typing import List
from collections import Counter, defaultdict

from utils_codon import is_valid_sequence

DP_MODULE = "DP"

class Decomposer:

    def __init__(self, mode=DP_MODULE):
        if mode not in ("DP",):
            # Only DP is supported in Codon-friendly port
            raise ValueError(f"{mode} is invalid mode for tandem repeat decomposer (Codon build supports DP only).")
        self.mode = "DP"

    @staticmethod
    def refine(decomposed_trs: List[List[str]], verbose: bool=False) -> List[List[str]]:
        """
            Resolve ambiguous boundaries by favoring the more frequent neighboring pair.
        """
        motif_pair_counter = Counter()
        motif_pair_str_counter = Counter()
        motif_pair_str_to_motif_pair = defaultdict(set)

        # Count motif pairs
        for tr in decomposed_trs:
            for i in range(len(tr) - 1):
                first_motif, second_motif = tr[i], tr[i+1]
                motif_pair = (first_motif, second_motif)
                motif_pair_str = ''.join(motif_pair)

                motif_pair_counter[motif_pair] += 1
                motif_pair_str_counter[motif_pair_str] += 1
                motif_pair_str_to_motif_pair[motif_pair_str].add(motif_pair)

        refined_trs: List[List[str]] = []
        for tr in decomposed_trs:
            new_tr: List[str] = []
            i = 0
            while i < len(tr):
                if i < len(tr) - 1:
                    first_motif, second_motif = tr[i], tr[i+1]
                    motif_pair = (first_motif, second_motif)

                    # After replacement, a new pair can be created. In this case, we just skip
                    if motif_pair_counter[motif_pair] == 0:
                        new_tr.append(first_motif)
                        i += 1
                        continue

                    motif_pair_str = ''.join(motif_pair)

                    best_pair = motif_pair
                    best_pair_count = motif_pair_counter[motif_pair]
                    for another_pair in motif_pair_str_to_motif_pair[motif_pair_str]:
                        another_pair_count = motif_pair_counter[another_pair]
                        if another_pair_count > best_pair_count:
                            best_pair = another_pair
                            best_pair_count = another_pair_count

                    new_tr.append(best_pair[0])
                    new_tr.append(best_pair[1])

                    # If we replaced, we need to decrement the counter and move to the next pair
                    motif_pair_counter[motif_pair] -= 1
                    motif_pair_counter[best_pair] += 1
                    i += 2
                else:
                    new_tr.append(tr[i])
                    i += 1

            refined_trs.append(new_tr)

        return refined_trs

    def decompose(self, sequence, motifs, **kwargs):
        """
        Decompose sequence into motifs using DP.
        """
        if not isinstance(sequence, str):
            raise TypeError("Sequence must be a string")
        if isinstance(motifs, str):
            motifs = [motifs]  # only one string

        sequence = sequence.upper()
        for i, motif in enumerate(motifs):
            motif = motif.upper()
            if not is_valid_sequence(motif):
                raise ValueError(f"Invalid character found in motif: {motif}")
            motifs[i] = motif

        if not is_valid_sequence(sequence):
            raise ValueError("Invalid character found in sequence")

        if self.mode == "DP":
            return self._decompose_dp(sequence, motifs, **kwargs)
        else:
            raise ValueError("Unsupported mode in Codon build")

    @staticmethod
    def _check_if_dp_parameters_are_valid(kwargs):
        # Provide defaults; accept overrides
        params = {
            "match_score": kwargs.get("match_score", 2.0),
            "mismatch_score": kwargs.get("mismatch_score", -1.0),
            "insertion_score": kwargs.get("insertion_score", -2.0),
            "min_score_threshold": kwargs.get("min_score_threshold", float("-inf")),
            "verbose": kwargs.get("verbose", False),
        }
        return params

    def _decompose_dp(self, sequence, motifs, **kwargs):
        params = self._check_if_dp_parameters_are_valid(kwargs)
        match_score = params["match_score"]
        mismatch_score = params["mismatch_score"]
        insertion_score = params["insertion_score"]
        min_score_threshold = params["min_score_threshold"]
        verbose = params["verbose"]

        n = len(sequence)
        M = len(motifs)
        Ls = [len(m) for m in motifs]
        maxL = max(Ls)

        # Create 3D tables using Python lists (Codon-friendly)
        # s[i][m][j]: best score for seq[:i] vs motif m up to j
        s = [[[float("-inf")] * (maxL + 1) for _ in range(M)] for _ in range(n + 1)]
        bt = {}

        # Initialize boundaries
        for m, motif in enumerate(motifs):
            s[0][m][0] = 0.0
            bt[(0,m,0)] = (0,m,0)
            # leading gaps along motif (deletions in motif)
            for j in range(1, len(motif) + 1):
                s[0][m][j] = s[0][m][j-1] + insertion_score
                bt[(0,m,j)] = (0,m,j-1)
            # leading gaps along sequence (insertions in motif)
            for i in range(1, n + 1):
                s[i][m][0] = s[i-1][m][0] + insertion_score
                bt[(i,m,0)] = (i-1,m,0)

        # Fill DP
        for i in range(1, n + 1):
            a = sequence[i-1]
            for m, motif in enumerate(motifs):
                L = len(motif)
                for j in range(1, L + 1):
                    b = motif[j-1]
                    # continue same motif: diag / up / left
                    from_diag = s[i-1][m][j-1] + (match_score if a == b else mismatch_score)
                    from_up   = s[i-1][m][j]   + insertion_score
                    from_left = s[i][m][j-1]   + insertion_score

                    best_score = from_diag
                    best_ptr = (i-1, m, j-1)
                    if from_up > best_score:
                        best_score, best_ptr = from_up, (i-1, m, j)
                    if from_left > best_score:
                        best_score, best_ptr = from_left, (i, m, j-1)

                    s[i][m][j] = best_score
                    bt[(i,m,j)] = best_ptr

                    # switching motifs when starting a motif (j==1)
                    if j == 1:
                        for pm, prev_motif in enumerate(motifs):
                            jump_src = (i-1, pm, len(prev_motif))
                            jump_val = s[jump_src[0]][jump_src[1]][jump_src[2]] + (match_score if a == b else mismatch_score)
                            if jump_val > s[i][m][j]:
                                s[i][m][j] = jump_val
                                bt[(i,m,j)] = jump_src

        # pick best end across motifs
        best_end = None
        best_val = float("-inf")
        for m, motif in enumerate(motifs):
            val = s[n][m][len(motif)]
            if val > best_val:
                best_val = val
                best_end = (n, m, len(motif))

        if best_val < min_score_threshold:
            if verbose:
                print("Best score below threshold:", best_val)
            return []

        # Backtrack to recover motif tokens
        path = []
        i, m, j = best_end
        while True:
            path.append((i,m,j))
            prev = bt[(i,m,j)]
            if (i,m,j) == prev:
                break
            i,m,j = prev
        path.reverse()

        tokens: List[str] = []
        curr_m = None
        curr_j_prev = None
        for (ii,mm,jj) in path[1:]:
            if curr_m is None:
                curr_m = mm
                curr_j_prev = jj
                continue
            if mm != curr_m:
                if curr_j_prev == len(motifs[curr_m]):
                    tokens.append(motifs[curr_m])
                curr_m = mm
                curr_j_prev = jj
            else:
                if jj == len(motifs[mm]) and curr_j_prev != jj:
                    tokens.append(motifs[mm])
                curr_j_prev = jj

        return tokens
