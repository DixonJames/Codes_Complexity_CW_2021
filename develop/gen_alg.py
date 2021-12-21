import numpy as np
import random
from ast import literal_eval
import time


population_size = 100
normal_mutation_const = 3
basis = None
best_tour = []




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
        try:
            norm = np.linalg.norm(components)
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




def mainSVP():
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


    coef_trial =[-2,  2, -1,  1, -3,  1 , 2 , 2, -2 , 1 , 2, -3]

    #72.0832851637604 [-2,  2, -1,  1, -3,  1 , 2 , 2, -2 , 1 , 2, -3] MSA

    shortest = LatticeVector(basis, coef_trial)
    vector = shortest.value()



class GenAlg:
    def __init__(self, time_frame, o_basis: Basis, original_basis: Basis, mutation_chance, pop_size, mutate=True):
        self.time_frame = time_frame
        self.o_basis = o_basis
        self.original_basis = original_basis
        self.mutation_chance = mutation_chance
        self.pop_size = pop_size
        self.mutate_option = mutate
    def max_tour(self, distances):
        total = 0
        for crossing in distances:
            total += max(crossing)
        return total

    def genStartPopulation(self, population_num, basis: Basis, minimum=0, maximum=1):
        population = []
        for p in range(population_num):
            rand_tour = basis.genRandPoint(minimum=minimum, maximum=maximum)
            population.append(rand_tour)
        return population

    def score(self, v: LatticeVector, original=True):
        #return v.norm()
        if original:
            s = v.changeBasis(self.original_basis, round=True).norm()
        else:
            s = v.norm()
        if s == 0.0:
            return 99999999999
        return s


    def findDuplicate(self, lst):
        dupes = []
        for char in lst:
            if lst.count(char) != 1:
                dupes.append(char)
        return dupes

    def basicCrossoverVectors(self, lattice, p, q):
        return lattice.modSivAvg(p, q, altBasis=self.original_basis), lattice.modSivDIff(p, q, altBasis=self.original_basis)

    def basicSwapMutation(self, v: LatticeVector):
        start = v.coeffs.copy()
        A_basis_vector = random.randint(0, len(v.coeffs) - 1)
        B_basis_vector = random.randint(0, len(v.coeffs) - 1)

        while B_basis_vector == A_basis_vector:
            B_basis_vector = random.randint(0, len(v.coeffs) - 1)

        v.coeffs[A_basis_vector], v.coeffs[B_basis_vector] = v.coeffs[B_basis_vector], v.coeffs[A_basis_vector]

        v.source = "swap"
        if v.nonZero():
            return v
        return self.basicSwapMutation(v)

    def tweakMutation(self, v:LatticeVector)-> LatticeVector:
        i = random.randint(0, len(v.coeffs) - 1)
        diff = random.choice([1, -1])

        new_v_coeffs = v.coeffs.copy()
        new_v_coeffs[i] = new_v_coeffs[i] + diff

        vn = LatticeVector(basis=v.basis, coeffs=new_v_coeffs)
        vn.source = "tweak"
        if vn.nonZero():
            return vn
        return self.tweakMutation(v)

    def applyMutations(self, new_pop, mutation_chance, mutate=True):
        final_new_pop = []
        if mutate==False:
            return new_pop
        for child in new_pop:
            if random.randint(1, mutation_chance) == 1:
                potencial_child = self.basicSwapMutation(child)

            elif random.randint(1, mutation_chance) == 1:
                potencial_child = self.tweakMutation(child)

            else:
                potencial_child = child


            final_new_pop.append(potencial_child)


        return final_new_pop

    def testPopulation(self, population, top_fitness, top_tour):
        for tour in population:
            if top_fitness >= self.score(tour):
                top_fitness = self.score(tour)
                top_tour = tour
        return top_fitness, top_tour

    def mixPops(self, old, new, basis):
        newlist = old + new + self.genStartPopulation(int(population_size/10), basis)
        old.extend(new)
        return sorted(newlist, key=self.score)[:population_size]

    def runTraining(self):
        start_time = time.time()
        tau = self.score(LatticeVector(self.o_basis, [3 for i in range(self.o_basis.dim)]))
        top_fitness = self.score(LatticeVector(self.o_basis, [10 for i in range(self.o_basis.dim)]))
        top_tour = []
        population = self.genStartPopulation(self.pop_size, self.o_basis, minimum=0, maximum=5)

        lattice = Lattice(self.o_basis)
        lattice.vectors = population

        elapsed_time = time.time() - start_time
        c_gen = 0

        tour_time_taken = time.time() - time.time()

        while (time.time() - start_time + tour_time_taken < self.time_frame):
            start_gen_time = time.time()

            population_fitness = [self.score(v) for v in lattice.vectors]
            total_fitness = sum(population_fitness)
            population_percentage = [1-((fitness)/ total_fitness)  for fitness in population_fitness]
            minFit = min(population_percentage)
            maxFit = max(population_percentage)

            population_percentage_a = np.array([1- (maxFit- fitness)/(maxFit-minFit) for fitness in population_percentage])
            population_percentage_b = np.array([1-((fitness) / max(population_fitness)) for fitness in population_fitness])
            population_percentage_c = np.array([((tau -   fitness) / tau) for fitness in population_fitness])
            controll =  np.array([1 for fitness in population_fitness])
            population_percentage = (population_percentage_b+ population_percentage_a + population_percentage_c)/3
            new_pop = []

            for j in range(self.pop_size):
                # making as many children as there are parents
                parents = random.choices(lattice.vectors, weights=population_percentage, k=2)
                #print(hash(parents[0]), hash(parents[1]))
                childA, childB = self.basicCrossoverVectors(lattice, parents[0], parents[1])

                childA_tour_length = self.score(childA)
                childB_tour_length = self.score(childB)

                if childA_tour_length <= childB_tour_length:
                    new_pop.append(childA)
                else:
                    new_pop.append(childB)

            lattice.vectors = self.mixPops(lattice.vectors, self.applyMutations(new_pop, self.mutation_chance), self.o_basis)

            top_fitness, top_tour = self.testPopulation(lattice.vectors, top_fitness, top_tour)
            print(top_tour.norm(), top_tour.coeffs, top_tour.source)
            c_gen += 1

            tour_time_taken = time.time() - start_gen_time

        # print(tour_time_taken, time.time() - start_time + tour_time_taken)
        return top_tour.norm(), top_tour.coeffs


def mainGA(normal_mutation_const=normal_mutation_const, population_size=population_size, mutate_option=False):
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
    path = "latticeBasis.txt"
    basis = Basis(basis_lst.copy())
    reduced_basis = LLL(Basis(basis_lst.copy()), delta=0.3).run()

    SVP_find = GenAlg(60, reduced_basis, basis, normal_mutation_const, population_size)
    tour, length = SVP_find.runTraining()
    return tour, length



if __name__ == '__main__':
    # naiveRun()
    for _ in range(100):
        print(mainGA(population_size = 60, mutate_option=False))