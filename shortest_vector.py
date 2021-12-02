import numpy as np
import math
import random
from ast import literal_eval


def genLines(file):
    while l := file.readline() != "":
        yield l


def readBasisFile(path):
    with open(path, 'r') as f:
        B = ""
        while (l := f.readline()) != "":
            if l[0] != "#" and l[0] != "B":
                B = B + l
    return np.array(literal_eval(B))


class LatticePoint:
    def __init__(self, basis, coeffs):
        self.basis = basis
        self.coeffs = coeffs

    def calcValue(self):
        """
        :return: co-oridnates of point
        """
        pass

    def __neg__(self):
        """
        :return: point with negative coeffcients
        """
        return LatticePoint(self.basis, [-(i) for i in self.coeffs])

    def __add__(self, other):
        return LatticePoint(self.basis, self.coeffs + other.coeffs)

class Basis:
    def __init__(self, vectors):
        self.vectors = vectors
        self.basisParts = self.vectors.T

    def genRandPoint(self):
        x = np.array([random.randint(0, len(self.vectors))])





def main():
    path = "latticeBasis.txt"
    basis = readBasisFile(path)


if __name__ == '__main__':
    main()
