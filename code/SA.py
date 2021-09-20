import sys
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio import AlignIO
from Bio.Align import substitution_matrices as sm
import numpy as np
import itertools


def align(seq1, seq2, **kwargs):
    match = kwargs.get("match", 2)
    mismatch = kwargs.get("mismatch", -1)
    gapOpening = kwargs.get("gapOpening", -0.5)
    gapExtending = kwargs.get("gapExtension", -0.1)
    alignments = pairwise2.align.globalms(seq1, seq2, match, mismatch, gapOpening, gapExtending)
    for a in [format_alignment(*a) for a in alignments]:
        print(a)


class Alignment:
    def __init__(self, file):
        temp = next(AlignIO.parse(file, "clustal"))
        self.records = temp._records
        self.alignment = np.array([r.seq for r in self.records])

    def getSeqID(self, ID):
        for r in self.records:
            if r.id == ID:
                return r.seq

    def getSeqPos(self, position):
        return self.records[position]

    def getCol(self, col):
        return self.alignment[:, col]

    def getPairScoreCol(self, colNum, mat, indel):
        column = self.getCol(colNum)
        combs = itertools.combinations(column, 2)
        matrix = sm.load(mat)
        score = 0
        for (x, y) in combs:
            if x == '-' and y == '-':
                continue
            elif x == '-' or y == '-':
                score += indel
            else:
                score += matrix.get((x, y))
        return score

    def getMSAScore(self, matrix, indel):
        score = 0
        aliLength = len(self.alignment[0])
        for i in range(0, aliLength):
            score += self.getPairScoreCol(i, matrix, indel)
        return score

    def getConservationScoreOnPos(self, pos):
        maximum = 0
        col = self.alignment.T[pos]
        for x in col:
            if x != '-':
                maximum = max(maximum, col.tolist().count(x))
        return maximum / len(col)

    def getTopNConserved(self, n):
        scores = []
        for x in range(0, self.alignment.T):
            scores.append(self.getConservationScoreOnPos(x))
        nBest = sorted(scores, reverse=True)[:n]
        conservedPositions = []
        i = -1
        for x in range(0, n):
            while i < len(scores):
                i += 1
                if scores[i] in nBest:
                    conservedPositions.append(i)
                    nBest.remove(scores[i])
                    break
        return conservedPositions
