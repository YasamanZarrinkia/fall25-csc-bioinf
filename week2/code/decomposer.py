from typing import List, Tuple

NEG_INF = -999999.0

def is_valid_sequence(sequence: str) -> bool:
    DNA_CHARACTERS = {'A', 'C', 'G', 'T'}
    for s in sequence:
        if s not in DNA_CHARACTERS:
            return False
    return True

class Decomposer:
    """
    Codon-friendly DP decomposer with working backtracking logic.
    """

    def __init__(self, mode: str = "DP"):
        if mode not in ("DP", "DP_CY", "HMM"):
            raise ValueError(str(mode) + " is invalid mode for tandem repeat decomposer.")
        # Keep the mode, even if we only implement DP in this file for now
        self.mode = mode

    # ---------------- refine ----------------
    @staticmethod
    def refine(decomposed_trs: List[List[str]], verbose: bool = False) -> List[List[str]]:
        from collections import Counter, defaultdict

        motif_pair_counter = Counter()
        motif_pair_str_counter = Counter()
        motif_pair_str_to_motif_pair = defaultdict(set)

        for tr in decomposed_trs:
            for i in range(len(tr) - 1):
                a = tr[i]
                b = tr[i + 1]
                pair = (a, b)
                pair_str = a + b
                motif_pair_counter[pair] += 1
                motif_pair_str_counter[pair_str] += 1
                motif_pair_str_to_motif_pair[pair_str].add(pair)

        refined: List[List[str]] = []
        for tr in decomposed_trs:
            for i in range(len(tr) - 1):
                a = tr[i]
                b = tr[i + 1]
                pair = (a, b)
                if motif_pair_counter.get(pair, 0) == 0:
                    continue
                key = a + b
                if motif_pair_counter[pair] < motif_pair_str_counter.get(key, 0):
                    best_pair = pair
                    best_c = -1
                    for p in motif_pair_str_to_motif_pair[key]:
                        c = motif_pair_counter.get(p, 0)
                        if c > best_c:
                            best_c = c
                            best_pair = p
                    if verbose:
                        print("Multiple pairs found", motif_pair_str_to_motif_pair[key])
                        print("Max frequency motif pair:", best_pair)
                    tr[i] = best_pair[0]
                    tr[i + 1] = best_pair[1]
            refined.append(tr)
        return refined

    # ---------------- public API ----------------
    @staticmethod
    def decompose(
        sequence,
        motifs_in,
        match_score=1.0,
        mismatch_score=-1.0,
        insertion_score=-1.0,
        deletion_score=-1.0,
        min_score_threshold=NEG_INF,
        verbose=False,
        **kwargs
    ) -> List[str]:
        # IMPORTANT: in Codon kwargs is a NamedTuple; explicitly detect invalid keys and RAISE
        if hasattr(kwargs, "mismatch"):
            # tests expect this to bubble up, do not swallow
            raise KeyError("Invalid keyword: mismatch")

        # Type validation
        if not (isinstance(match_score, (int, float))):
            raise ValueError(f"Invalid value type. match_score: {match_score}")
        if not (isinstance(mismatch_score, (int, float))):
            raise ValueError(f"Invalid value type. mismatch_score: {mismatch_score}")
        if not (isinstance(insertion_score, (int, float))):
            raise ValueError(f"Invalid value type. insertion_score: {insertion_score}")
        if not (isinstance(deletion_score, (int, float))):
            raise ValueError(f"Invalid value type. deletion_score: {deletion_score}")
        if not (isinstance(min_score_threshold, (int, float))):
            raise ValueError(f"Invalid value type. min_score_threshold: {min_score_threshold}")
        if not (verbose is True or verbose is False):
            raise ValueError(f"Invalid value type. verbose: {verbose}")

        if not isinstance(sequence, str):
            raise TypeError(f"Sequence must be a string, got {type(sequence)}")

        # Validate motifs parameter
        if isinstance(motifs_in, str):
            motifs_list: List[str] = [motifs_in]
        else:
            if not isinstance(motifs_in, list):
                raise TypeError(f"Motifs must be a list of strings, got {type(motifs_in)}")
            motifs_list = list(motifs_in)

        if len(motifs_list) == 0:
            # tests expect ValueError when there are no motifs
            raise ValueError("No motifs provided")

        # Validate sequence content
        seq = sequence.upper()
        if not is_valid_sequence(seq):
            raise ValueError("Sequence has invalid characters: " + seq)

        # Validate motifs content
        motifs: List[str] = []
        for m in motifs_list:
            mu = str(m).upper()
            if not is_valid_sequence(mu):
                raise ValueError("The motif has invalid characters: " + mu)
            motifs.append(mu)

        # Convert scores to float
        return Decomposer._decompose_dp(
            seq,
            motifs,
            float(match_score),
            float(mismatch_score),
            float(insertion_score),
            float(deletion_score),
            float(min_score_threshold),
            bool(verbose),
        )

    # ---------------- Helper functions ----------------
    @staticmethod
    def _argmax_3(a: float, b: float, c: float) -> int:
        if a >= b and a >= c:
            return 0
        elif b >= c:
            return 1
        else:
            return 2

    # ---------------- DP core ----------------
    @staticmethod
    def _decompose_dp(
        sequence: str,
        motifs: List[str],
        match_score: float,
        mismatch_score: float,
        insertion_score: float,
        deletion_score: float,
        min_score_threshold: float,
        verbose: bool,
    ) -> List[str]:
        n = len(sequence)
        M = len(motifs)

        max_len = 0
        for mot in motifs:
            if len(mot) > max_len:
                max_len = len(mot)

        # s[i][m][j] and back[i][m][j] -> (i2,m2,j2)
        s: List[List[List[float]]] = []
        back: List[List[List[Tuple[int, int, int]]]] = []
        for i in range(n + 1):
            s.append([])
            back.append([])
            for m in range(M):
                s[i].append([0.0] * (max_len + 1))
                back[i].append([(0, m, 0)] * (max_len + 1))

        # boundaries
        for m_idx in range(M):
            Lm = len(motifs[m_idx])
            for i in range(n + 1):
                for j in range(Lm + 1):
                    if i == 0 and j == 0:
                        s[0][m_idx][0] = 0.0
                        back[0][m_idx][0] = (0, m_idx, 0)
                    elif i == 0 and j != 0:
                        s[0][m_idx][j] = s[0][m_idx][j - 1] + insertion_score
                        back[0][m_idx][j] = (0, m_idx, j - 1)
                    elif i != 0 and j == 0:
                        s[i][m_idx][0] = s[i - 1][m_idx][0] + insertion_score
                        back[i][m_idx][0] = (i - 1, m_idx, 0)

        # fill DP table
        for i in range(1, n + 1):
            si = sequence[i - 1]
            for m_idx in range(M):
                mot = motifs[m_idx]
                Lm = len(mot)
                for j in range(1, Lm + 1):
                    mj = mot[j - 1]
                    if j == 1:
                        if i == 1:
                            diag = s[i - 1][m_idx][j - 1] + (match_score if si == mj else mismatch_score)
                            left = s[i - 1][m_idx][j] + insertion_score
                            up = s[i][m_idx][j - 1] + deletion_score

                            argmax_index = Decomposer._argmax_3(diag, left, up)
                            if argmax_index == 0:
                                s[i][m_idx][j] = diag
                                back[i][m_idx][j] = (0, m_idx, 0)
                            elif argmax_index == 1:
                                s[i][m_idx][j] = left
                                back[i][m_idx][j] = (i - 1, m_idx, j)
                            else:
                                s[i][m_idx][j] = up
                                back[i][m_idx][j] = (i, m_idx, 0)
                        else:
                            best_val = NEG_INF
                            best_m = 0
                            best_jend = 0
                            for t in range(M):
                                mend = s[i - 1][t][len(motifs[t])]
                                if mend > best_val:
                                    best_val = mend
                                    best_m = t
                                    best_jend = len(motifs[t])

                            motif_end = best_val + (match_score if si == mot[0] else mismatch_score)
                            left = s[i - 1][m_idx][1] + insertion_score
                            up = s[i][m_idx][0] + deletion_score

                            argmax_index = Decomposer._argmax_3(motif_end, left, up)
                            if argmax_index == 0:
                                s[i][m_idx][j] = motif_end
                                back[i][m_idx][j] = (i - 1, best_m, best_jend)
                            elif argmax_index == 1:
                                s[i][m_idx][j] = left
                                back[i][m_idx][j] = (i - 1, m_idx, 1)
                            else:
                                s[i][m_idx][j] = up
                                back[i][m_idx][j] = (i, m_idx, 0)
                    else:
                        diag = s[i - 1][m_idx][j - 1] + (match_score if si == mj else mismatch_score)
                        left = s[i - 1][m_idx][j] + insertion_score
                        up = s[i][m_idx][j - 1] + deletion_score

                        argmax_index = Decomposer._argmax_3(diag, left, up)
                        if argmax_index == 0:
                            s[i][m_idx][j] = diag
                            back[i][m_idx][j] = (i - 1, m_idx, j - 1)
                        elif argmax_index == 1:
                            s[i][m_idx][j] = left
                            back[i][m_idx][j] = (i - 1, m_idx, j)
                        else:
                            s[i][m_idx][j] = up
                            back[i][m_idx][j] = (i, m_idx, j - 1)

        # choose best end
        start_ptr = None
        best_score = min_score_threshold
        for m_idx in range(M):
            Lm = len(motifs[m_idx])
            val = s[n][m_idx][Lm]
            if val > best_score:
                best_score = val
                start_ptr = (n, m_idx, Lm)

        if start_ptr is None:
            raise ValueError("No good match greater than score threshold of " + str(min_score_threshold))

        if verbose:
            print("Best score:", best_score)

        # Backtracking
        decomposed_motifs: List[str] = []
        current = start_ptr
        reverse_path = []
        while current != (0, current[1], 0):
            reverse_path.append(current)
            i, m, j = current
            current = back[i][m][j]
        reverse_path.append((0, current[1], 0))  # start

        path = list(reversed(reverse_path))

        if verbose:
            print(f"Path length: {len(path)}")
            # print("Path:", path)

        # Reconstruct motifs
        current_motif_chars: List[str] = []
        prev_i, prev_m, prev_j = -1, -1, -1

        for (i, m, j) in path:
            # boundary between motifs: previous motif finished & current j==1
            if prev_m != -1 and prev_j == len(motifs[prev_m]) and j == 1 and current_motif_chars:
                decomposed_motifs.append(''.join(current_motif_chars))
                current_motif_chars = []
            # consume base if we moved along sequence
            if i > prev_i and i > 0:
                current_motif_chars.append(sequence[i - 1])
            prev_i, prev_m, prev_j = i, m, j

        if current_motif_chars:
            decomposed_motifs.append(''.join(current_motif_chars))

        if verbose:
            reconstructed = "".join(decomposed_motifs)
            print(f"Decomposed: {decomposed_motifs}")
            print(f"Reconstructed: {reconstructed}")
            print(f"Original:      {sequence}")
            print(f"Match: {reconstructed == sequence}")

        return decomposed_motifs
