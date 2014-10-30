# requires scikit bio [http://scikit-bio.org]

from skbio.parse.sequences import parse_fasta
import sys

if __name__ == '__main__':
    with open(sys.argv[1], 'U') as in_file:
        with open(sys.argv[2], 'w') as out_file:
            for label, seq in parse_fasta(in_file):
                seq = seq.replace("U", "T")
                out_file.write(">%s\n%s\n" % (label, seq))
