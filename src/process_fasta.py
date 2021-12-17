"""Process PASTA files."""
import sys
import numpy as np


def parse_fasta(
        path: str,
) -> (list, list):
    """Parse headers and sequences of a PASTA file.

    Args:
        path: path to PASTA file

    Returns:
        tuple of list of headers and list of sequences
    """
    content: list = open(path, 'r').readlines()
    headers: list = list()
    sequences: list = list()
    current_seq: str = ""

    for line in content:
        if line.startswith('>'):  # new header
            headers.append((line[1:]).strip())
            if len(current_seq) > 0:  # current sequence not empty
                sequences.append(current_seq)
                current_seq = ""
        else:  # append line of continuous sequence
            current_seq += line.strip()

    # account for last sequence
    if len(current_seq) > 0:
        sequences.append(current_seq)

    return headers, sequences


def discard_ambiguous_seqs(
        content: list,
) -> list:
    """Discard any sequences with non DNA base letters.

    Args:
        content: list of sequences

    Returns:
        list of sequences which exclusively consist of DNA base letters
    """
    unambiguous_seqs: list = list()
    seq_index: int = 0

    while seq_index < len(content):
        seq_uppercase: str = content[seq_index].upper()
        valid_counter: int = 0

        for c in seq_uppercase:
            if c not in ['A', 'C', 'G', 'T']:
                break
            valid_counter += 1

        if valid_counter == len(content[seq_index]):
            unambiguous_seqs.append(content[seq_index])

        seq_index += 1

    return unambiguous_seqs


def nucleotide_frequencies(
        content: list,
) -> None:
    """Print frequencies of nucleotides.

    Args:
        content: list of sequences

    Returns:
        None
    """
    counters: np.array = np.array([0, 0, 0, 0])  # A,C,G,T

    for seq in content:
        for c in seq:
            if c in ['a', 'A']:
                counters[0] += 1
            elif c in ['c', 'C']:
                counters[1] += 1
            elif c in ['g', 'G']:
                counters[2] += 1
            elif c in ['t', 'T']:
                counters[3] += 1

    total_nucs: np.array = np.sum(counters)
    frequencies: np.array = np.round(counters / total_nucs, 1)

    print('A:', frequencies[0])
    print('C:', frequencies[1])
    print('G:', frequencies[2])
    print('T:', frequencies[3])


def map_reads(
        query_path: str,
        reference_path: str
) -> dict:
    """Map reads to reference sequences.

    Parse FASTA files, print nucleotide frequencies, and map query sequences to
    given reference sequences.

    Args:
        query_path: path to FASTA file with query sequences
        reference_path: path to FASTA file with reference sequences

    Returns:
        dictionary of query sequences and dictionary of reference sequences and
        list of matching positions in reference sequence
    """
    # parse files
    ambiguous_query_names, ambiguous_query_seqs = parse_fasta(query_path)
    reference_names, reference_seqs = parse_fasta(reference_path)

    # discard queries with non-DNA letters
    query_seqs: list = discard_ambiguous_seqs(ambiguous_query_seqs)
    query_names: list = list()
    i: int = 0
    for seq in range(len(ambiguous_query_seqs)):
        if ambiguous_query_seqs[seq].__eq__(query_seqs[i]):
            query_names.append(ambiguous_query_names[seq])
            i += 1

    # print nucleotide fractions
    print('Nucleotide fractions queries:')
    nucleotide_frequencies(query_seqs)
    print('Nucleotide fractions references:')
    nucleotide_frequencies(reference_seqs)

    # create dictionaries
    outer_dict: dict = {}
    pos: int

    for q in range(len(query_seqs)):
        inner_dict = {}

        for r in range(len(reference_seqs)):
            start_positions = list()
            pos = reference_seqs[r].upper().find(query_seqs[q].upper()) + 1

            while pos > 0:
                start_positions.append(pos)
                pos = reference_seqs[r].upper().find(query_seqs[q].upper(),
                                                     pos) + 1

            # if len(start_positions) > 0:
            inner_dict[reference_names[r]] = start_positions

        outer_dict[query_names[q]] = inner_dict

    return outer_dict


if __name__ == '__main__':
    try:
        map_reads(sys.argv[1], sys.argv[2])
    except IndexError or FileNotFoundError:
        print('ERROR: Please input "python process_fasta.py '
              '<path to query file> <path to reference file>"!')
