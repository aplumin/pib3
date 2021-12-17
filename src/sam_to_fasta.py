"""
Convert SAM to FASTA file.
"""
import sys


def main(sam_path, fasta_path):
    """Convert a SAM file to a FASTA file.

    Args:
        sam_path: path to SAM file
        fasta_path: path to FASTA file

    Returns:
        None
    """
    content = open(sam_path, 'r').readlines()
    headers = list()
    sequences = list()

    for line in content:
        line_fields = line.split('\t')
        current_header = ""
        if not line.startswith('@'):
            if len(line) > 0:
                for field in range(9):
                    current_header += line_fields[field]
                    current_header += " "
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
    print("SUCCESS: FASTA file", sys.argv[2], "created!")


if __name__ == '__main__':
    try:
        main(sys.argv[1], sys.argv[2])
    except:
        print('ERROR: Please input "python sam_to_fasta.py <path to existing SAM file> <path to FASTA file>"!')
