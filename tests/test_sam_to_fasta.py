"""Test sam_to_fasta."""
import pytest
from _pytest.capture import CaptureFixture
from src.sam_to_fasta import main


def create_sam_valid(
        path: str,
) -> None:
    """Create a SAM file for testing.

    Args:
        path: path to test SAM file

    Returns:
        None
    """
    test_file = open(path, 'w')
    test_file.write("@Header 42 \n@more\theader\tinfo\n")
    test_file.write("QNAME\tFLAG\tRNAME\tPOS\tMAPQ\tCIGAR"
                    "\tRNEXT\tPNEXT\tTLEN\t")
    test_file.write("SOME DNA SEQUENCE\tQUAL\n")
    test_file.write("qname\tflag\trname\tpos\tmapq\tcigar"
                    "\trnext\tpnext\ttlen\t")
    test_file.write("sequence\tqual\toptional field"
                    "\tother optional field\n")
    test_file.close()


def create_sam_invalid(
        path: str,
) -> None:
    """Create a file not conforming to SAM formatting.

    Args:
        path: path to invalid SAM file

    Returns:
        None
    """
    test_file = open(path, 'w')
    test_file.write("SOME CONTENT")
    test_file.close()


def test_sam_to_fasta_valid(
) -> None:
    """Test whether a valid SAM file is correctly converted to a FASTA file.

    Returns:
        None
    """
    create_sam_valid("test_valid.sam")
    main("test_valid.sam", "test_valid.fa")
    fasta_content = open("test_valid.fa", 'r').readlines()
    assert len(fasta_content) == 4
    assert fasta_content[0] == ">QNAME\tFLAG\tRNAME\tPOS\tMAPQ\tCIGAR" \
                               "\tRNEXT\tPNEXT\tTLEN\n"
    assert fasta_content[1] == "SOME DNA SEQUENCE\n"
    assert fasta_content[2] == ">qname\tflag\trname\tpos\tmapq\tcigar" \
                               "\trnext\tpnext\ttlen\n"
    assert fasta_content[3] == "sequence\n"


def test_sam_to_fasta_invalid_path(
) -> None:
    """Test whether FileNotFoundError if an invalid path is given.

    Returns:
        None
    """
    # Invalid
    with pytest.raises(FileNotFoundError):
        main("invalid_path.sam", "test_valid.fa")


def test_sam_to_fasta_invalid_format(
        capfd: CaptureFixture,
) -> None:
    """Test whether an incorrect format produces correct errors.

    Args:
        capfd: CaptureFixture of output

    Returns:
        None
    """
    create_sam_invalid("test_invalid.sam")
    main("test_invalid.sam", "test_invalid.fa")
    # FASTA not created
    with pytest.raises(FileNotFoundError):
        open("test_invalid.fa", 'r')
    # Output
    out, err = capfd.readouterr()
    assert out == "ERROR: not a valid SAM format!\n"
