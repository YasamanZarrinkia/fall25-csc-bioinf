import string
from typing import Dict, List, Optional, Tuple

LOWERCASE_LETTERS = string.ascii_lowercase
UPPERCASE_LETTERS = string.ascii_uppercase
DIGITS = string.digits

skipping_characters = ['(', '=', '<', '>', '?', '-']
PRIVATE_MOTIF_LABEL = '?'
INDEX_TO_CHR = list(LOWERCASE_LETTERS) + list(UPPERCASE_LETTERS) + list(DIGITS)
INDEX_TO_CHR.extend([chr(x) for x in range(33, 127) if chr(x) not in skipping_characters and chr(x) not in INDEX_TO_CHR])

def get_motif_counter(decomposed_vntrs: List[List[str]]) -> Dict[str, int]:
    """Return a counter for each motif using plain dict"""
    motif_counter: Dict[str, int] = {}
    for decomposed_vntr in decomposed_vntrs:
        for motif in decomposed_vntr:
            if motif in motif_counter:
                motif_counter[motif] += 1
            else:
                motif_counter[motif] = 1
    return motif_counter

class MotifEncoder:
    """
    Codon-friendly MotifEncoder:
      - no self-imports
      - top-level class definition
      - concrete dict types (no Any/object)
    """

    def __init__(self, private_motif_threshold: int = 0):
        self.private_motif_threshold: int = private_motif_threshold
        self.symbol_to_motif: Dict[str, str] = {}
        self.motif_to_symbol: Dict[str, str] = {}
        self.score_matrix: Dict[str, Dict[str, float]] = {}
        self.motif_counter: Dict[str, int] = {}

    @staticmethod
    def _divide_motifs_into_normal_and_private(
        motif_counter: Dict[str, int],
        private_motif_threshold: int
    ) -> Tuple[Dict[str, int], Dict[str, int]]:
        normal: Dict[str, int] = {}
        private: Dict[str, int] = {}
        items = list(motif_counter.items())
        items.sort(key=lambda kv: (-kv[1], kv[0]))
        for motif, cnt in items:
            if cnt > private_motif_threshold:
                normal[motif] = cnt
            else:
                private[motif] = cnt
        return normal, private

    @staticmethod
    def find_private_motif_threshold(
        decomposed_vntrs: List[List[str]],
        label_count: Optional[int] = None
    ) -> int:
        max_labels = len(INDEX_TO_CHR) - 1
        if label_count is not None:
            max_labels = label_count - 1  # 1 reserved for private
        motif_counter = get_motif_counter(decomposed_vntrs)
        items = list(motif_counter.items())
        items.sort(key=lambda kv: (-kv[1], kv[0]))
        threshold = 0
        used = 0
        for _motif, cnt in items:
            used += 1
            if used > max_labels:
                threshold = cnt
                break
        return threshold

    @staticmethod
    def write_motif_map(
        output_file: str,
        motif_to_alphabet: Dict[str, str],
        motif_counter: Dict[str, int]
    ) -> None:
        items = list(motif_counter.items())
        items.sort(key=lambda kv: (-kv[1], kv[0]))
        with open(output_file, "w") as f:
            for motif, _ in items:
                f.write(motif + "\t" + str(motif_to_alphabet[motif]) + "\t" + str(motif_counter[motif]) + "\n")

    @staticmethod
    def _encode_decomposed_tr(
        decomposed_vntrs: List[List[str]],
        motif_to_symbol: Dict[str, str]
    ) -> List[str]:
        out: List[str] = []
        for vntr in decomposed_vntrs:
            parts: List[str] = []
            for m in vntr:
                parts.append(motif_to_symbol[m])
            out.append("".join(parts))
        return out

    def encode(
        self,
        decomposed_vntrs: List[List[str]],
        motif_map_file: str,
        label_count: Optional[int] = None,
        auto: bool = True
    ) -> List[str]:
        def idx_to_char(i: int) -> str:
            if i < 0 or i >= len(INDEX_TO_CHR):
                raise ValueError("Index should range 0.." + str(len(INDEX_TO_CHR)-1) + ", got " + str(i))
            return INDEX_TO_CHR[i]

        # compute threshold
        if label_count is not None:
            self.private_motif_threshold = self.find_private_motif_threshold(decomposed_vntrs, label_count)
        if auto:
            self.private_motif_threshold = self.find_private_motif_threshold(decomposed_vntrs)

        # motif counts
        motif_counter = get_motif_counter(decomposed_vntrs)
        self.motif_counter = {}
        for k, v in motif_counter.items():
            self.motif_counter[k] = v

        motif_to_symbol: Dict[str, str] = {}
        symbol_to_motif: Dict[str, str] = {}

        if self.private_motif_threshold > 0:
            normal, private = self._divide_motifs_into_normal_and_private(motif_counter, self.private_motif_threshold)
            if len(normal) + 1 > len(INDEX_TO_CHR):
                raise ValueError("Too many unique motifs: " + str(len(normal) + len(private)))

            # assign shared label to all private motifs
            for motif in private:
                motif_to_symbol[motif] = PRIVATE_MOTIF_LABEL
            # representative back-map (last one wins)
            symbol_to_motif[PRIVATE_MOTIF_LABEL] = ""
            for motif in private:
                symbol_to_motif[PRIVATE_MOTIF_LABEL] = motif

            # assign labels to normal motifs by frequency
            items = list(normal.items())
            items.sort(key=lambda kv: (-kv[1], kv[0]))
            i = 0
            for motif, _cnt in items:
                ch = idx_to_char(i)
                motif_to_symbol[motif] = ch
                symbol_to_motif[ch] = motif
                i += 1
        else:
            if len(motif_counter) > len(INDEX_TO_CHR):
                raise ValueError("Too many unique motifs: " + str(len(motif_counter)))
            items = list(motif_counter.items())
            items.sort(key=lambda kv: (-kv[1], kv[0]))
            i = 0
            for motif, _cnt in items:
                ch = idx_to_char(i)
                motif_to_symbol[motif] = ch
                symbol_to_motif[ch] = motif
                i += 1

        # persist map
        self.write_motif_map(motif_map_file, motif_to_symbol, motif_counter)

        # store maps
        self.motif_to_symbol = {}
        for k, v in motif_to_symbol.items():
            self.motif_to_symbol[k] = v
        self.symbol_to_motif = {}
        for k, v in symbol_to_motif.items():
            self.symbol_to_motif[k] = v

        # (score_matrix not used in current Codon pipeline; leave empty dict)
        self.score_matrix = {}

        return self._encode_decomposed_tr(decomposed_vntrs, motif_to_symbol)
