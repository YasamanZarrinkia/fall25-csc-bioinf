from collections import Counter
import string

LOWERCASE_LETTERS = string.ascii_lowercase
UPPERCASE_LETTERS = string.ascii_uppercase
DIGITS = string.digits

skipping_characters = ['(', '=', '<', '>', '?', '-']
PRIVATE_MOTIF_LABEL = '?'
INDEX_TO_CHR = list(LOWERCASE_LETTERS) + list(UPPERCASE_LETTERS) + list(DIGITS)
INDEX_TO_CHR.extend([chr(x) for x in range(33, 127) if chr(x) not in skipping_characters and chr(x) not in INDEX_TO_CHR])

DNA_CHARACTERS = {'A', 'C', 'G', 'T'}

def get_motif_counter(decomposed_vntrs: list[list[str]]) -> Counter:
    """Return a counter for each motif"""
    motif_counter = Counter()
    for decomposed_vntr in decomposed_vntrs:
        motif_counter.update(Counter(decomposed_vntr))
    return motif_counter

def is_valid_sequence(sequence: str) -> bool:
    """Check if the given sequence is DNA sequence"""
    for s in sequence:
        if s not in DNA_CHARACTERS:
            return False
    return True

def get_levenshtein_distance(s1: str, s2: str) -> int:
    """
    Calculate Levenshtein distance between two strings
    """
    if len(s1) < len(s2):
        return get_levenshtein_distance(s2, s1)

    if len(s2) == 0:
        return len(s1)

    previous_row = list(range(len(s2) + 1))
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous_row[j + 1] + 1
            deletions = current_row[j] + 1
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row

    return previous_row[-1]

def get_score_matrix(symbol_to_motif: dict,
                     match_score: int = 2,
                     mismatch_score_for_edit_dist_of_1: int = -1,
                     mismatch_score_for_edit_dist_greater_than_1: int = -2,
                     gap_open_penalty: float = 1.5,
                     gap_extension_penalty: float = 0.6) -> dict:
    """
    Create score matrix for alignment
    """
    score_matrix = {}
    score_matrix['gap_open'] = gap_open_penalty
    score_matrix['gap_extension'] = gap_extension_penalty

    for symbol1 in symbol_to_motif:
        score_matrix[symbol1] = {}
        for symbol2 in symbol_to_motif:
            motif_seq1 = symbol_to_motif[symbol1]
            motif_seq2 = symbol_to_motif[symbol2]
            if symbol1 == symbol2:
                score_matrix[symbol1][symbol2] = match_score
            else:
                edit_dist = get_levenshtein_distance(motif_seq1, motif_seq2)
                edit_dist_cutoff = 1
                if abs(len(motif_seq1) - len(motif_seq2)) <= 1:
                    edit_dist_cutoff += len(max(motif_seq2, motif_seq1, key=len)) // 30
                if edit_dist <= edit_dist_cutoff:
                    score_matrix[symbol1][symbol2] = mismatch_score_for_edit_dist_of_1
                else:
                    score_matrix[symbol1][symbol2] = mismatch_score_for_edit_dist_greater_than_1

    return score_matrix

def add_padding(encoded_trs: list[str]) -> list[str]:
    """
    Add padding to encoded TRs to make them equal length
    """
    max_motif_count = len(max(encoded_trs, key=len))
    padded_trs = []
    for encoded_tr in encoded_trs:
        padding_count = max_motif_count - len(encoded_tr)
        padded_trs.append(encoded_tr + '-' * padding_count)
    return padded_trs

def sort(aligned_vntrs: list[str], sample_ids: list[str], symbol_to_motif: dict, 
         sample_order_file: str = None, method: str = 'motif_count') -> tuple[list[str], list[str]]:
    """Sort the aligned and encoded tandem repeats"""
    if method == 'name':
        sorted_data = sorted(zip(sample_ids, aligned_vntrs), key=lambda x: x[0])
        return [x[0] for x in sorted_data], [x[1] for x in sorted_data]
    elif method == 'motif_count':
        sorted_data = sorted(zip(sample_ids, aligned_vntrs), 
                           key=lambda x: len(x[1].replace('-', '')), reverse=True)
        return [x[0] for x in sorted_data], [x[1] for x in sorted_data]
    elif method == 'manually':
        if sample_order_file is None:
            return sample_ids, aligned_vntrs
        # Simplified manual sorting - in practice you'd read from file
        return sample_ids, aligned_vntrs
    else:
        raise ValueError(f"Please check the rearrangement method. {method}")
