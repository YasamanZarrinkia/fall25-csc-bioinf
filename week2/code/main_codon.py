
import sys
from typing import List
from decomposer_codon import Decomposer
from motif_encoder_codon import MotifEncoder
from motif_aligner_codon import MotifAligner
from utils_codon import sort, add_padding, get_motif_marks

class TandemRepeatVizWorker:
    def __init__(self):
        self.decomposer = Decomposer()
        self.motif_encoder = MotifEncoder()
        self.motif_aligner = MotifAligner()

    def generate_trplot(self,
                        tr_id: str,
                        sample_ids: List[str],
                        tr_sequences: List[str],
                        motifs: List[str],
                        skip_alignment: bool=False,
                        rearrangement_method: str='name',
                        sample_order_file: str=None,
                        output_dir: str='./',
                        pad_left: int = 0,
                        pad_right: int = 0,
                        show_figure: bool=False,
                        output_name: str=None,
                        **kwargs):
        # 1) decompose
        decomposed = [self.decomposer.decompose(seq, motifs, **kwargs) for seq in tr_sequences]
        # 2) refine
        decomposed = self.decomposer.refine(decomposed)
        # 3) encode
        encoded_vntrs, symbol_to_motif, score_matrix, motif_counter = self.motif_encoder.encode(decomposed, score_matrix=None)
        # 4) align
        if not skip_alignment:
            sample_ids, encoded_vntrs = self.motif_aligner.align(sample_ids, encoded_vntrs, tr_id, score_matrix, output_dir, tool='mafft')
        # padding
        encoded_vntrs = add_padding(encoded_vntrs, pad_left, pad_right)
        # 5) sort
        sample_ids, encoded_vntrs = sort(sample_ids, encoded_vntrs, method=rearrangement_method, sample_order_file=sample_order_file)
        # NOTE: Visualization is intentionally omitted in Codon build. Use Python visualizer on saved outputs.
        return sample_ids, encoded_vntrs, symbol_to_motif, score_matrix, motif_counter

if __name__ == "__main__":
    # Tiny CLI for smoke test:
    # python main_codon.py ATG,TTG S1:S2:S3 ATGATGTTG:ATGATGATG:ATGTTGTTG
    if len(sys.argv) != 4:
        print("Usage: main_codon.py <comma_motifs> <colon_sample_ids> <colon_sequences>")
        sys.exit(1)
    motifs = sys.argv[1].split(',')
    sample_ids = sys.argv[2].split(':')
    sequences = sys.argv[3].split(':')
    w = TandemRepeatVizWorker()
    sids, aln, sym2motif, score, counts = w.generate_trplot("demo", sample_ids, sequences, motifs)
    print("SAMPLES:", sids)
    print("ALN:")
    for sid, row in zip(sids, aln):
        print(sid, row)
    print("SYMBOLS→MOTIFS:", sym2motif)
    print("COUNTS:", dict(counts))

