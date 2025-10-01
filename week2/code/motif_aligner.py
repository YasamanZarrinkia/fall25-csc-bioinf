class MotifAligner:
    
    def align(self, sample_ids: list[str], encoded_vntrs: list[str], vid: str = None,
              score_matrix: dict = None, output_dir: str = "./", tool: str = "center_star") -> tuple[list[str], list[str]]:
        """
        Align encoded VNTRs using center-star alignment
        """
        if tool == "center_star":
            return self._align_motifs_with_center_star(sample_ids, encoded_vntrs)
        else:
            raise ValueError(f"Unsupported alignment tool: {tool}")

    @staticmethod
    def _align_motifs_with_center_star(sample_ids: list[str], encoded_vntrs: list[str]) -> tuple[list[str], list[str]]:
        """
        Perform center-star alignment - simplified version
        """
        if len(encoded_vntrs) == 0:
            return sample_ids, encoded_vntrs

        # Use first sequence as center
        center_seq = encoded_vntrs[0]
        center_id = sample_ids[0]
        
        aligned_seqs = [center_seq]
        aligned_ids = [center_id]

        # Align all other sequences to center with simple gap insertion
        for i in range(1, len(encoded_vntrs)):
            seq = encoded_vntrs[i]
            sample_id = sample_ids[i]
            
            # Simple alignment: if sequences have different lengths, pad with gaps
            max_len = max(len(center_seq), len(seq))
            padded_center = center_seq + '-' * (max_len - len(center_seq))
            padded_seq = seq + '-' * (max_len - len(seq))
            
            aligned_seqs.append(padded_seq)
            aligned_ids.append(sample_id)

        return aligned_ids, aligned_seqs

    @staticmethod
    def load_aligned_trs(aln_output: str) -> tuple[list[str], list[str]]:
        """
        Load aligned TRs from file (simplified - in practice you'd parse alignment file)
        """
        # This is a placeholder - in real implementation you'd parse the alignment file
        aligned_trs = []
        aligned_sample_ids = []
        
        try:
            with open(aln_output, "r") as handle:
                for line in handle:
                    if line.startswith(">"):
                        aligned_sample_ids.append(line[1:].strip())
                    else:
                        aligned_trs.append(line.strip())
        except:
            # Fallback if file doesn't exist
            pass
            
        return aligned_sample_ids, aligned_trs
