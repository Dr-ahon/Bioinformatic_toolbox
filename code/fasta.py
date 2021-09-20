import sys
import argparse

from Bio import SeqIO


def get_generator(file):
    return SeqIO.parse(file, "fasta").records


def desc(name, file):
    for record in file:
        if record.name == name:
            return record.description


def length(name, file):
    for record in file:
        if record.name == name:
            return len(record.seq)


def sequence(name, file):
    for record in file:
        if record.name == name:
            return record.seq


def subsequence(name, file, start, end):
    for record in file:
        if record.name == name:
            return record.seq[start: end]


def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", action="store", dest="fileName")
    parser.add_argument("-d-", "--desc", action="store", dest="descName")
    parser.add_argument("-q", "--seq", action="store", dest="seqName")
    parser.add_argument("-l", "--length", action="store", dest="lenName")
    parser.add_argument("-s", "--subseq", nargs=3, action="store", dest="subseqArgs")

    args = parser.parse_args()
    if args.fileName:
        fileName = args.fileName
    else:
        print("filename missing")
        sys.exit()

    with open(fileName) as file:

        if args.descName:
            print(desc(args.descName, SeqIO.parse(file, "fasta")))

        elif args.seqName:
            print(sequence(args.seqName, SeqIO.parse(file, "fasta")))
            file = SeqIO.parse(file, "fasta")
        elif args.lenName:
            print(length(args.lenName, SeqIO.parse(file, "fasta")))

        elif args.subseqArgs:
            print(subsequence(args.subseqArgs[0], SeqIO.parse(file, "fasta"), int(args.subseqArgs[1]),
                              int(args.subseqArgs[2])))


if __name__ == '__main__':
    main(sys.argv[1:])
