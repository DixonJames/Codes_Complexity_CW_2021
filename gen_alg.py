import random
import time
from shortest_vector import *

population_size = 100
normal_mutation_const = 3
basis = None
best_tour = []


class GenAlg:
    def __init__(self, time_frame, o_basis: Basis, original_basis: Basis, mutation_chance, pop_size):
        self.time_frame = time_frame
        self.o_basis = o_basis
        self.original_basis = original_basis
        self.mutation_chance = mutation_chance
        self.pop_size = pop_size

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
        return v.norm()
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

    def applyMutations(self, new_pop, mutation_chance):
        final_new_pop = []
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
            #print(top_tour.norm(), top_tour.coeffs, top_tour.source)
            c_gen += 1

            tour_time_taken = time.time() - start_gen_time

        # print(tour_time_taken, time.time() - start_time + tour_time_taken)
        return top_tour.norm(), top_tour.coeffs


def main(normal_mutation_const=normal_mutation_const, population_size=population_size):
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

    SVP_find = GenAlg(60, basis, basis, normal_mutation_const, population_size)
    tour, length = SVP_find.runTraining()
    return tour, length

if __name__ == '__main__':
    # naiveRun()
    for population in range(50, 100, 10):
      for i in range(5):
        print(population,  main(population_size = population))