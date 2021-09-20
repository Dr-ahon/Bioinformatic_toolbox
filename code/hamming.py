import argparse
import sys
from Bio import SeqIO
import fasta


def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", action="store", dest="fileName")
    parser.add_argument("-s", "--seqs", action="store", dest="seqNames", nargs=2)

    args = parser.parse_args()

    if args.fileName:
        if args.seqNames:
            seqName1 = args.seqNames[0]
            seqName2 = args.seqNames[1]
            seq1 = fasta.sequence(seqName1, SeqIO.parse(args.fileName, "fasta"))
            seq2 = fasta.sequence(seqName2, SeqIO.parse(args.fileName, "fasta"))

            if len(seq1) == len(seq2):
                count = 0
                for i, j in zip(seq1, seq2):
                    if i != j: count += 1
                print(count)
            else:
                print("sequences are of different lengths")
        else:
            print("missing sequence")
    else:
        print("missing fileName")


if __name__ == '__main__':
    main(sys.argv[1:])
