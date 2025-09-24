
from collections import Counter
from typing import Dict, List

from utils_codon import INDEX_TO_CHR, PRIVATE_MOTIF_LABEL
from utils_codon import get_motif_counter, get_score_matrix

class MotifEncoder:

    def __init__(self, private_motif_threshold=0):
        self.private_motif_threshold = private_motif_threshold
        self.symbol_to_motif = None
        self.motif_to_symbol = None
        self.score_matrix = None
        self.motif_counter = None

    @staticmethod
    def _divide_motifs_into_normal_and_private(motif_counter, private_motif_threshold):
        normal_motifs = Counter()
        private_motifs = Counter()
        for motif, cnt in motif_counter.items():
            if cnt > private_motif_threshold:
                normal_motifs[motif] = cnt
            else:
                private_motifs[motif] = cnt
        return normal_motifs, private_motifs

    def encode(self, decomposed_vntrs: List[List[str]], score_matrix=None):
        motif_counter = get_motif_counter(decomposed_vntrs)
        normal_motifs, private_motifs = self._divide_motifs_into_normal_and_private(motif_counter, self.private_motif_threshold)

        # assign symbols to normal motifs
        symbol_to_motif = dict()
        motif_to_symbol = dict()
        i = 0
        for motif, _ in normal_motifs.items():
            symbol = INDEX_TO_CHR[i]
            symbol_to_motif[symbol] = motif
            motif_to_symbol[motif] = symbol
            i += 1

        # private motifs map to '?'
        for motif, _ in private_motifs.items():
            motif_to_symbol[motif] = PRIVATE_MOTIF_LABEL

        # encode
        encoded_vntrs = []
        for vntr in decomposed_vntrs:
            encoded_vntrs.append(''.join([motif_to_symbol[m] for m in vntr]))

        # score matrix
        if score_matrix is None:
            score_matrix = get_score_matrix(symbol_to_motif)

        self.symbol_to_motif = symbol_to_motif
        self.motif_to_symbol = motif_to_symbol
        self.score_matrix = score_matrix
        self.motif_counter = motif_counter

        return encoded_vntrs, symbol_to_motif, score_matrix, motif_counter

