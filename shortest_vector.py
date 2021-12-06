import numpy as np
import math
import random
from ast import literal_eval


def readBasisFile(path):
    with open(path, 'r') as f:
        B = ""
        while (l := f.readline()) != "":
            if l[0] != "#" and l[0] != "B":
                B = B + l
    return np.array(literal_eval(B))


def randPairing(num):
    nums = [i for i in range(num)]
    pairs = []
    for i in range(num // 2):
        ca = random.choice(nums)
        nums.remove(ca)
        cb = random.choice(nums)
        nums.remove(cb)
        pairs.append((ca, cb))
    return pairs


class LatticeVector:
    def __init__(self, basis, coeffs):
        self.basis = basis
        self.coeffs = coeffs

    def value(self):
        """
        :return: co-oridnates of point
        """
        try:
            return np.multiply(self.coeffs, self.basis.vectors)
        except:
            print("d")

    def norm(self):
        components = self.value()
        t = 0
        sqr = lambda x: x ** 2
        norm = math.sqrt(sum(sqr(sum(components.T))))
        return norm

    def __neg__(self):
        """
        :return: point with negative coeffcients
        """
        return LatticeVector(self.basis, [-(i) for i in self.coeffs])

    def __add__(self, other):
        return LatticeVector(self.basis, self.coeffs + other.coeffs)

    def __sub__(self, other):
        s = LatticeVector(self.basis, [-(i) for i in other.coeffs])
        return LatticeVector(self.basis, self.coeffs + s.coeffs)

    def innerProd(self, other) -> int:
        return sum(self.basis * other.basis)


class Basis:
    def __init__(self, vectors):
        self.vectors = vectors
        self.dim = len(self.vectors)

        self.orthogonal = None

    def genRandPoint(self):
        x = np.array(
            [random.randint(0, math.ceil(self.dim + math.log(self.dim, 2))) for _ in range(self.dim)])
        return LatticeVector(self, x)

    def createOrthoginal(self):
        pass


class GramSchmidt:
    def __init__(self, original_basis):
        self.original_basis = original_basis
        self.orthogonal_basis = self.run()

    def proj(self, u, v, normalised=False):
        u,v = np.array(u), np.array(v)
        if not normalised:
            return sum(u*v.T) * u
        return (sum(u*v.T)/sum(u*u.T)) * u.coeffs

    def run(self):
        us = [list(self.original_basis.vectors[0])]

        for k in range(1, len(self.original_basis.vectors)):
            u_k = self.original_basis.vectors[k]

            for j in range(k - 1):
                u_k = u_k - self.proj(us[j], self.original_basis.vectors[k])

            us.append(u_k)

        return np.array(us)


class Lattice:
    def __init__(self, basis):
        self.basis = basis
        self.vectors = []

        self.norms = dict()

    def addVector(self, v: LatticeVector):
        check = sum(v.coeffs % 1)
        if check == 0:
            if f"{v.coeffs}" not in self.norms.keys():
                self.vectors.append(v)
                self.norms[f"{v.coeffs}"] = v.norm()

    def rankVectors(self):
        s_lst = sorted(self.vectors, key=lambda v: v.norm())
        return s_lst

    def sivDiff(self, va: LatticeVector, vb: LatticeVector) -> LatticeVector:
        """
        sieving by differences
        :param va: LatticeVector
        :param vb: LatticeVector
        :return: LatticeVector
        """
        return va - vb

    def sivAvg(self, va: LatticeVector, vb: LatticeVector) -> LatticeVector:
        """
        sieving by averages
        :param va:
        :param vb:
        :return:
        """
        coeffs = (va.coeffs + vb.coeffs) / 2
        check = sum(coeffs % 1)
        if check != 0:
            return LatticeVector(self, coeffs)
        return None

    def modSivAvg(self, va: LatticeVector, vb: LatticeVector) -> LatticeVector:
        """
        modified sieving by differences
        :param va:
        :param vb:
        :return:
        """
        coeffs = (va.coeffs + vb.coeffs) / 2

        f = lambda x: x * random.choice([2, -2])
        modification = [f(x) for x in coeffs % 1]

        n_vb = LatticeVector(self.basis, vb.coeffs + modification)
        coeffs = (va.coeffs + n_vb.coeffs) / 2
        coeffs = coeffs.astype(int)

        return LatticeVector(self.basis, coeffs)


class SVP:
    def __init__(self, basis):
        self.basis = basis
        self.lattice = Lattice(self.basis)

    def initialise(self, number_points):
        for _ in range(number_points):
            v = self.basis.genRandPoint()
            self.lattice.addVector(v)

    def naiveCombo(self, const_num):
        self.initialise(const_num)

        for i in range(100):
            r = self.lattice.rankVectors()
            top = r
            pairs = randPairing(const_num)
            new = [self.lattice.modSivAvg(top[pair[0]], top[pair[1]]) for pair in pairs]

            self.lattice.vectors = r[int(const_num * 0.1):] + new

        return self.lattice.rankVectors()[0]


def testDiffs():
    path = "latticeBasis.txt"
    basis = Basis(readBasisFile(path))
    lattice = Lattice(basis)

    p = basis.genRandPoint()
    q = basis.genRandPoint()
    lattice.vectors.append(p)
    lattice.vectors.append(q)

    diff = lattice.sivDiff(p, q)
    avg = lattice.sivAvg(p, q)
    m_avg = lattice.modSivAvg(p, q)

    n = p.norm()
    print(n)

def naiveRun():
    path = "latticeBasis.txt"
    basis = Basis(readBasisFile(path))
    svp = SVP(basis)

    shortest_vec = svp.naiveCombo(1000)
    print(shortest_vec.norm())

def main():
    path = "latticeBasis.txt"
    basis = Basis(readBasisFile(path))
    orth = GramSchmidt(basis).original_basis


if __name__ == '__main__':
    main()
