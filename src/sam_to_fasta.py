"""Convert SAM to FASTA file."""
import sys


def main(
        sam_path: str,
        fasta_path: str,
) -> None:
    """Convert a SAM file to a FASTA file.

    Args:
        sam_path: path to SAM file
        fasta_path: path to FASTA file

    Returns:
        None
    """
    content: list = open(sam_path, 'r').readlines()
    headers: list = list()
    sequences: list = list()

    for line in content:
        line_fields: list = line.split('\t')
        current_header: str = ""
        if not line.startswith('@'):
            if len(line) > 0:
                try:
                    for field in range(8):
                        current_header += line_fields[field]
                        current_header += "\t"
                except IndexError:
                    print("ERROR: not a valid SAM format!")
                    return
                current_header += line_fields[8]
                headers.append(current_header)
                sequences.append(line_fields[9])

    fasta = open(fasta_path, 'w')
    for seq in range(len(headers)):
        fasta.write(">")
        fasta.write(headers[seq])
        fasta.write("\n")
        fasta.write(sequences[seq])
        fasta.write("\n")
    fasta.close()
    print("SUCCESS: FASTA file", fasta_path, "created!")


if __name__ == '__main__':
    try:
        main(sys.argv[1], sys.argv[2])
    except IndexError or FileNotFoundError:
        print('ERROR: Please input as arguments '
              '"<path to existing SAM file> <path to FASTA file>"!')
