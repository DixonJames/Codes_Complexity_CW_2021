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


class Basis:
    def __init__(self, vectors):
        self.vectors = vectors
        self.dim = len(self.vectors)

        self.orthogonal = None

    def genRandPoint(self, minimum=0, maximum=1):
        x = np.array([0 for _ in range(self.dim)])
        while not any(x != np.array([0 for _ in range(self.dim)])):
            x = np.array([random.randint(minimum, maximum) for _ in range(self.dim)])
        v = LatticeVector(self, x)
        v.source = "rand"
        return v

    def createOrthoginal(self):
        self.orthogonal = Basis(GramSchmidt(self).orthogonal_basis)


class LatticeVector:
    def __init__(self, basis, coeffs):
        self.basis = basis
        self.coeffs = coeffs

        self.coeff_identity = np.array([0 for _ in range(len(coeffs))])
        self.source = None

    def __hash__(self):
        return int(str(hash(str(self.coeffs)))[:5])

    def nonZero(self, altBasis=None):
        if altBasis is None:
            if not all(self.coeffs == self.coeff_identity):
                return True
            return False

        alotbasis_v = self.changeBasis(altBasis)
        if not all(self.coeffs == self.coeff_identity) and not all(alotbasis_v.coeffs == alotbasis_v.coeff_identity):
            return True
        return False

    def value(self):
        """
        :return: co-oridnates of point
        """
        try:
            return np.array([sum(v) for v in np.multiply(self.basis.vectors, self.coeffs)])
        except:
            print("d")

    def norm(self):
        components = self.value()
        t = 0
        sqr = lambda x: x ** 2
        try:
            norm = np.linalg.norm(components.T)
        except:
            print("s")
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

    def changeBasis(self, new_basis: Basis, round=True):
        """
        :param new_basis: basis to convert into
        :return: coefficents  for self's vector for the input basis
        """
        U = self.basis.vectors
        W = new_basis.vectors
        inv_U = np.linalg.inv(U)
        inv_W = np.linalg.inv(W)

        if round:
            nb_coeffs = (np.matmul(inv_W, np.matmul(U, self.coeffs))).astype(int)
            b_coeffs = (np.matmul(inv_U, np.matmul(W, nb_coeffs))).astype(int)
        else:
            nb_coeffs = (np.matmul(inv_W, np.matmul(U, self.coeffs)))
            b_coeffs = (np.matmul(inv_U, np.matmul(W, nb_coeffs)))

        return LatticeVector(new_basis, coeffs=nb_coeffs)


class GramSchmidt:
    def __init__(self, original_basis):
        self.original_basis = original_basis
        self.orthogonal_basis = self.run()

    def proj(self, u, v, normalised=False):
        u, v = np.array(u), np.array(v)
        if not normalised:
            return np.dot(u, v) * u

        return np.dot(u, v) / np.dot(u, u) * u

    def run(self):
        us = np.array([list(self.original_basis.vectors[0])])

        for k in range(1, len(self.original_basis.vectors)):
            u_k = self.original_basis.vectors[k]

            basis_component = u_k - sum(self.proj(u, self.original_basis.vectors[v]) for u, v in
                                        zip([u for u in us], [k for _ in range(k)]))

            us = np.append(us, [basis_component], axis=0)

        return us


class LLL:
    def __init__(self, basis: Basis, delta=0.26):
        self.basis = basis
        self.basis.createOrthoginal()
        self.delta = delta

    def dot(self, u, v):
        return np.dot(u.T, v.T)

    def createU(self, i, j):
        basis_dim = self.basis.dim
        return self.dot(self.basis.vectors[i], self.basis.orthogonal.vectors[j]) / self.dot(
            self.basis.orthogonal.vectors[j], self.basis.orthogonal.vectors[j])

    def run(self):
        """
        based on pseudocode fromT theorem 6.68 in:
        @book{silverman2008introduction,
        title={An introduction to mathematical cryptography},
        author={Silverman, Joseph H and Pipher, Jill and Hoffstein, Jeffrey},
        year={2008},
        publisher={Springer}
        }
        :return: a reduced basis
        """
        k = 1
        while k < self.basis.dim:
            for j in range(k - 1, -1, -1):
                if abs(self.createU(k, j)) > 0.5:
                    self.basis.vectors[k] = self.basis.vectors[k] - (self.basis.vectors[j] * int(self.createU(k, j)))
                    self.basis.createOrthoginal()

            k_self_dot = self.dot(self.basis.orthogonal.vectors[k], self.basis.orthogonal.vectors[k])
            kBelow_self_dot = self.dot(self.basis.orthogonal.vectors[k - 1], self.basis.orthogonal.vectors[k - 1])
            if k_self_dot >= (self.delta - self.createU(k, k - 1) ** 2) * kBelow_self_dot:
                k += 1
            else:
                self.basis.vectors[k], self.basis.vectors[k - 1] = self.basis.vectors[k - 1], self.basis.vectors[k]
                self.basis.createOrthoginal()
                k = max(1, k - 1)
        return self.basis


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
        res = va - vb
        if res.nonZero():
            return res
        return va

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
            res = LatticeVector(self, coeffs)
        if res.nonZero():
            return res
        return va

    def modSivAvg(self, va: LatticeVector, vb: LatticeVector, altBasis=None) -> LatticeVector:
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

        res = LatticeVector(self.basis, coeffs)

        res.source = "MSA"
        if res.nonZero(altBasis):
            return res
        return va

    def modSivDIff(self, va: LatticeVector, vb: LatticeVector, altBasis=None):
        diff = self.sivDiff(va, vb)
        s = diff
        for candidate_vec in self.vectors:
            c_diff = self.sivDiff(diff, candidate_vec)
            if c_diff.norm() < diff.norm() and c_diff.nonZero(altBasis):
                s = c_diff
        s.source = "MSD"
        return s


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
            pairs = randPairing(min(const_num, len(r)))
            new = [self.lattice.modSivAvg(top[pair[0]], top[pair[1]]) for pair in pairs]

            self.lattice.vectors = r[int(const_num * 0.1):] + new
            print(self.lattice.rankVectors()[0].norm())

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


def basisChangeTest():
    path = "latticeBasis.txt"
    basis = Basis(readBasisFile(path))
    basis.createOrthoginal()
    orth = basis.orthogonal

    va = LatticeVector(basis, [1 for i in range(basis.dim)])
    vb = va.changeBasis(orth, round=False)

    vc = vb.changeBasis(basis, round=False)


def main():
    path = "latticeBasis.txt"
    basis_lst = np.array([[37, 20, 96, 20, 34, 64, 82, 56, 47, 21, 50, 49],
                          [39, 24, 19, 49, 82, 97, 88, 84, 41, 51, 36, 74],
                          [19, 56, 37, 73, 4, 12, 72, 18, 46, 8, 54, 94],
                          [13, 46, 26, 8, 83, 71, 45, 84, 21, 32, 53, 80],
                          [65, 39, 25, 56, 52, 44, 84, 30, 69, 33, 13, 5],
                          [59, 56, 90, 1, 42, 58, 90, 92, 2, 6, 7, 80],
                          [18, 14, 26, 31, 91, 93, 77, 64, 95, 36, 23, 5],
                          [11, 58, 22, 51, 90, 13, 93, 43, 21, 81, 12, 77],
                          [42, 65, 99, 6, 23, 43, 94, 30, 37, 66, 34, 66],
                          [99, 31, 24, 44, 18, 58, 17, 27, 70, 88, 59, 11],
                          [30, 43, 21, 70, 48, 47, 13, 93, 94, 48, 69, 58],
                          [7, 12, 94, 88, 59, 95, 43, 62, 71, 36, 91, 70]])
    basis = Basis(basis_lst)
    basis.createOrthoginal()
    reduces_basis = LLL(Basis(readBasisFile(path))).run()

    coef_trial = [0, 0, 0, -1, 1, -1, 0, 0, 0, 0, 1, 0]  # 97.0
    coef_trial = [0, -1, 0, 0, 0, -1, 1, 1, -1, 0, 2, -1]  # 80.8
    # 76.87652437513027 [ 2 -2  0 -1  2 -1 -1 -2  2 -1 -1  3] MSA
    # 74.96665925596525 [ 1 -1  0 -1  1  0 -1 -1  1  0 -1  2] MSD
    #72.0832851637604 [-2,  2, -1,  1, -3,  1 , 2 , 2, -2 , 1 , 2, -3] MSA

    va = LatticeVector(basis, coef_trial)
    vb = va.changeBasis(reduces_basis, round=False)


if __name__ == '__main__':
    # naiveRun()
    main()
