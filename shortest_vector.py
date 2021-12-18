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

    def genRandPoint(self):
        x = np.array([0 for _ in range(self.dim)])
        while not any(x != np.array([0 for _ in range(self.dim)])):
            x = np.array([random.randint(0, 3) for _ in range(self.dim)])
        return LatticeVector(self, x)

    def createOrthoginal(self):
        self.orthogonal = Basis(GramSchmidt(self).orthogonal_basis)


class LatticeVector:
    def __init__(self, basis, coeffs):
        self.basis = basis
        self.coeffs = coeffs

        self.coeff_identity = np.array([0 for _ in range(len(coeffs))])

    def nonZero(self):
        if not all(self.coeffs == self.coeff_identity):
            return True
        return False

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

    def changeBasis(self, new_basis: Basis, round=True):
        """
        :param new_basis: basis to convert into
        :return: coefficents  for self's vector for the input basis
        """
        U = self.basis.vectors.T
        W = new_basis.vectors.T
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
    def __init__(self, basis: Basis, delta=3 / 4):
        self.orginal_basis = basis
        self.delta = delta

        self.orginal_basis.createOrthoginal()
        self.original_orthogonal_basis = self.orginal_basis.orthogonal

        self.u = self.createU(self.original_orthogonal_basis, self.orginal_basis)

    def createU(self, orthogonal, original):
        basis_dim = self.orginal_basis.dim
        return [[np.dot(b_i, b_j) for b_j in orthogonal] for b_i in original.vectors]

    def run(self):
        b_copy = self.orginal_basis
        o_b_copy = self.original_orthogonal_basis
        u_copy = self.u

        k = 1
        while k <= self.orginal_basis.dim:
            for j in range(k - 1, 0, -1):
                if u_copy[k][j] > 1 / 2:
                    b_copy.vectors[k] = o_b_copy[k] - u_copy[k][j] * b_copy.vectors[j]

                    # the update
                    b_copy.createOrthoginal()
                    o_b_copy = self.orginal_basis.orthogonal
                    u_copy = self.createU(o_b_copy, b_copy)

            if np.dot(o_b_copy[k], o_b_copy[k]) >= (self.delta - (u_copy[k][k - 1]) ** 2) * np.dot(o_b_copy[k - 1],
                                                                                                   o_b_copy[k - 1]):
                k += 1
            else:
                b_copy.vectors[k], b_copy.vectors[k - 1] = b_copy.vectors[k - 1], b_copy.vectors[k]
                # the update
                b_copy.createOrthoginal()
                o_b_copy = self.orginal_basis.orthogonal
                u_copy = self.createU(o_b_copy, b_copy)
                k = max(k - 1, 1)

        return b_copy


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

        res = LatticeVector(self.basis, coeffs)

        if res.nonZero():
            return res
        return va

    def modSivDIff(self, va: LatticeVector, vb: LatticeVector):
        diff = self.sivDiff(va, vb)
        s = diff
        for candidate_vec in self.vectors:
            c_diff = self.sivDiff(diff, candidate_vec)
            if c_diff.norm() < diff.norm():
                s = c_diff
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
    basis.createOrthoginal()
    orth = basis.orthogonal

    va = LatticeVector(basis, [1 for i in range(basis.dim)])
    vb = va.changeBasis(orth, round=False)

    vc = vb.changeBasis(basis, round=False)

    reduces_basis = LLL(basis).run()




if __name__ == '__main__':
    main()
