"""Test process_fasta."""

from pathlib import Path

import pytest
from _pytest.capture import CaptureFixture

from src.process_fasta import parse_fasta, discard_ambiguous_seqs,\
    nucleotide_frequencies, map_reads


def create_fasta(
        path: str,
) -> (list, list):
    """Create a test FASTA file.

    Args:
        path: path to test FASTA file

    Returns:
        tuple of FASTA list with headers list with contents
    """
    # Create test file
    test_file = open(path, 'w')
    test_file.write(">seq1  \nTEST\n>seq2\nSEQUENCE\n  OV ER"
                    "\nSEVERAL \nLINES")
    test_file.close()
    return parse_fasta(path)


def test_parse_fasta_length(
) -> None:
    """Test whether the header and content lists have correct lengths.

    Returns:
        None
    """
    h, s = create_fasta(
        Path(__file__).parents[0] / 'test_files' / 'test_parse_fasta.fa',
    )
    assert len(h) == 2
    assert len(s) == 2


def test_parse_fasta_content(
) -> None:
    """Test content of header and content lists.

    Returns:
        None
    """
    h, s = create_fasta(
        Path(__file__).parents[0] / 'test_files' / 'test_parse_fasta.fa',
    )
    assert h[0] == "seq1"
    assert h[1] == "seq2"
    assert s[0] == "TEST"
    assert s[1] == "SEQUENCEOVERSEVERALLINES"


def test_parse_fasta_invalid(
) -> None:
    """Test whether FileNotFoundError is raised with an invalid path.

    Returns:
        None
    """
    with pytest.raises(FileNotFoundError):
        parse_fasta("a")


def test_discard_ambiguous_seqs_empty(
) -> None:
    """Test whether empty sequences are not discarded.

    Returns:
        None
    """
    assert discard_ambiguous_seqs([]) == []
    assert discard_ambiguous_seqs(["", ""]) == ["", ""]


def test_discard_ambiguous_seqs_valid(
) -> None:
    """Test whether valid DNA sequences are not discarded.

    Returns:
        None
    """
    assert discard_ambiguous_seqs(["", "AGCCTT", "agcctt", "AgCcTT"]) \
           == ["", "AGCCTT", "agcctt", "AgCcTT"]


def test_discard_ambiguous_seqs_invalid(
) -> None:
    """Test whether invalid DNA sequences are discarded.

    Returns:
        None
    """
    assert discard_ambiguous_seqs([" ", "agt", "42", "car", "@?0\n"]) \
           == ["agt"]
    with pytest.raises(TypeError):
        discard_ambiguous_seqs(42)


def test_nucleotide_frequencies_empty(
        capfd: CaptureFixture,
) -> None:
    """Test whether sequences without nucleotide letters are identified.

    Args:
        capfd: CaptureFixture of output

    Returns:
        None
    """
    nucleotide_frequencies([])
    out, err = capfd.readouterr()
    assert out == "A: 0\nC: 0\nG: 0\nT: 0\n"

    nucleotide_frequencies(["number 42!"])
    out, err = capfd.readouterr()
    assert out == "A: 0\nC: 0\nG: 0\nT: 0\n"


def test_nucleotide_frequencies_valid(
        capfd: CaptureFixture,
) -> None:
    """Test whether frequencies with nucleotide letters are printed correctly.

    Args:
        capfd: CaptureFixture of output

    Returns:
        None
    """
    nucleotide_frequencies(["A"])
    out, err = capfd.readouterr()
    assert out == "A: 1.0\nC: 0.0\nG: 0.0\nT: 0.0\n"

    nucleotide_frequencies(["AAAAACCCGGT"])
    out, err = capfd.readouterr()
    assert out == "A: 0.45\nC: 0.27\nG: 0.18\nT: 0.09\n"


def test_nucleotide_frequencies_mixed(
        capfd: CaptureFixture,
) -> None:
    """Test frequencies of sequences also containing other letters.

    Args:
        capfd: CaptureFixture of output

    Returns:
        None
    """
    nucleotide_frequencies(["AAXAAACCCGGTX"])
    out, err = capfd.readouterr()
    assert out == "A: 0.45\nC: 0.27\nG: 0.18\nT: 0.09\n"


def init_map_reads(
        query_path_1: str,
        query_path_2: str,
        reference_path: str,
) -> None:
    """Create test files for test_map_reads.

    Args:
        query_path_1: path to test query file 1
        query_path_2: path to test query file 2
        reference_path: path to test reference file

    Returns:
        None
    """
    query_file = open(query_path_1, 'w')
    query_file.write(">seq1\nAAAAA\n>seq2\nCCCCC")
    query_file.close()
    query_file_2 = open(query_path_2, 'w')
    query_file_2.write(">seq1\nAA\n>seq2\nAG")
    query_file_2.close()
    reference_file = open(reference_path, 'w')
    reference_file.write(">ref1\nAAAAA\n>ref2\nTTTTT")
    reference_file.close()


def test_map_reads_output(
        capfd: CaptureFixture,
) -> None:
    """Test whether query sequences are mapped correctly to references.

    Args:
        capfd: CaptureFixture of output

    Returns:
        None
    """
    init_map_reads(
        Path(__file__).parents[0] / 'test_files' / 'test_query_1.fa',
        Path(__file__).parents[0] / 'test_files' / 'test_query_2.fa',
        Path(__file__).parents[0] / 'test_files' / 'test_reference.fa',
    )
    map_reads(
        Path(__file__).parents[0] / 'test_files' / 'test_query_1.fa',
        Path(__file__).parents[0] / 'test_files' / 'test_reference.fa',
    )
    out, err = capfd.readouterr()
    assert out == "Nucleotide fractions queries:\n" \
                  "A: 0.5\nC: 0.5\nG: 0.0\nT: 0.0\n" \
                  "Nucleotide fractions references:\n" \
                  "A: 0.5\nC: 0.0\nG: 0.0\nT: 0.5\n"


def test_map_reads_dictionary(
) -> None:
    """Test dictionary that map_reads returns.

    Returns:
        None
    """
    init_map_reads(
        Path(__file__).parents[0] / 'test_files' / 'test_query_1.fa',
        Path(__file__).parents[0] / 'test_files' / 'test_query_2.fa',
        Path(__file__).parents[0] / 'test_files' / 'test_reference.fa',
    )
    read_dict_1 = map_reads(
        Path(__file__).parents[0] / 'test_files' / 'test_query_1.fa',
        Path(__file__).parents[0] / 'test_files' / 'test_reference.fa',
    )
    assert len(read_dict_1) == 2
    assert read_dict_1["seq1"]["ref1"] == [1]  # offset 1
    assert read_dict_1["seq1"]["ref2"] == []
    assert read_dict_1["seq2"]["ref1"] == []
    assert read_dict_1["seq2"]["ref2"] == []

    # multiple hits in one sequence
    read_dict_2 = map_reads(
        Path(__file__).parents[0] / 'test_files' / 'test_query_2.fa',
        Path(__file__).parents[0] / 'test_files' / 'test_reference.fa',
    )
    assert len(read_dict_2) == 2
    assert read_dict_2["seq1"]["ref1"] == [1, 2, 3, 4]  # offset 1
    assert read_dict_2["seq1"]["ref2"] == []
    assert read_dict_2["seq2"]["ref1"] == []
    assert read_dict_2["seq2"]["ref2"] == []
