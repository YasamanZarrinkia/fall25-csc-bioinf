from decomposer import Decomposer
from motif_encoder import MotifEncoder
from motif_aligner import MotifAligner

class TandemRepeatVizWorker:

    def __init__(self):
        self.decomposer = Decomposer()
        self.motif_encoder = MotifEncoder()
        self.motif_aligner = MotifAligner()

    def generate_trplot(self,
                        tr_id: str,
                        sample_ids: list,
                        tr_sequences: list,
                        motifs: list,
                        skip_alignment: bool = False,
                        rearrangement_method: str = 'motif_count',
                        sample_order_file: str = None,
                        output_dir: str = "./",
                        verbose: bool = True):
        """
        Main workflow: decomposition → encoding → alignment → rearrangement
        """
        if len(sample_ids) != len(tr_sequences):
            raise ValueError("The number of sample IDs and sequences are different.")

        if verbose:
            print(f"ID: {tr_id}")
            print(f"Motifs: {motifs}")
            print(f"Loaded {len(tr_sequences)} tandem repeat sequences")
            print("Decomposing TR sequences")

        # 1. Decomposition
        decomposed_trs = []
        for i, tr_sequence in enumerate(tr_sequences):
            if verbose:
                print(f"Progress: {i+1}/{len(tr_sequences)}")
            decomposed_trs.append(self.decomposer.decompose(tr_sequence, motifs))

        # Refinement of decomposition
        decomposed_trs = self.decomposer.refine(decomposed_trs)

        # 2. Encoding
        if verbose:
            print("Encoding")
        encoded_trs = self.motif_encoder.encode(decomposed_trs,
                                              motif_map_file=f"{output_dir}/{tr_id}_motif_map.txt",
                                              auto=True)

        # 3. Alignment
        if skip_alignment:
            if verbose:
                print("Skip alignment step")
            # Simple padding
            max_len = max(len(seq) for seq in encoded_trs)
            aligned_trs = [seq + '-' * (max_len - len(seq)) for seq in encoded_trs]
            sorted_sample_ids = sample_ids
        else:
            if verbose:
                print("Alignment using center-star")
            sorted_sample_ids, aligned_trs = self.motif_aligner.align(sample_ids, encoded_trs, tr_id,
                                                                     output_dir=output_dir, tool='center_star')

        if verbose:
            print("Processing completed")

        return sorted_sample_ids, aligned_trs, decomposed_trs
