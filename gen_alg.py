import random
import time
from shortest_vector import *

population_size = 50
normal_mutation_const = 2

best_tour = []


def max_tour(distances):
    total = 0
    for crossing in distances:
        total += max(crossing)
    return total


def genStartPopulation(population_num, basis: Basis):
    population = []
    for p in range(population_num):
        rand_tour = basis.genRandPoint()
        population.append(rand_tour)
    return population


def score(v: LatticeVector):
    return v.norm()



def findDuplicate(lst):
    dupes = []
    for char in lst:
        if lst.count(char) != 1:
            dupes.append(char)
    return dupes


def basicCrossoverVectors(lattice, p, q):
    return lattice.modSivAvg(p, q), lattice.modSivDIff(p, q)


def basicSwapMutation(v: LatticeVector):
    A_basis_vector = random.randint(0, len(v.coeffs) - 1)
    B_basis_vector = random.randint(0, len(v.coeffs) - 1)

    while B_basis_vector == A_basis_vector:
        B_basis_vector = random.randint(0, len(v.coeffs) - 1)

    v.coeffs[A_basis_vector], v.coeffs[B_basis_vector] = v.coeffs[B_basis_vector], v.coeffs[A_basis_vector]

    if v.nonZero():
        return v
    return None


def applyMutations(new_pop, mutation_chance):
    final_new_pop = []
    for child in new_pop:
        if random.randint(1, mutation_chance) == 1:
            new_child = basicSwapMutation(child)

        else:
            new_child = child
        final_new_pop.append(new_child)

    return final_new_pop


def testPopulation(population, top_fitness, top_tour):
    for tour in population:
        if top_fitness >= score(tour):
            top_fitness = score(tour)
            top_tour = tour
    return top_fitness, top_tour


def mixPops(old, new, basis):
    newlist = old + new + genStartPopulation(10, basis)
    old.extend(new)
    return sorted(old, key=score)[:100]


def runTraining(time_frame, basis: Basis, mutation_chance, pop_size):
    start_time = time.time()
    tau = LatticeVector(basis, [3 for i in range(basis.dim)]).norm()
    top_fitness = LatticeVector(basis, [3 for i in range(basis.dim)]).norm()
    top_tour = []
    population = genStartPopulation(pop_size, basis)

    lattice = Lattice(basis)
    lattice.vectors = population

    elapsed_time = time.time() - start_time
    c_gen = 0

    tour_time_taken = time.time() - time.time()

    while (time.time() - start_time + tour_time_taken < time_frame):
        start_gen_time = time.time()

        population_fitness = [score(v) for v in lattice.vectors]
        total_fitness = sum(population_fitness)
        population_percentage = [(tau - fitness) / total_fitness for fitness in population_fitness]

        new_pop = []

        for j in range(pop_size):
            # making as many children as there are parents
            parents = random.choices(lattice.vectors, weights=population_percentage, k=2)
            childA, childB = basicCrossoverVectors(lattice, parents[0], parents[1])

            childA_tour_length = score(childA)
            childB_tour_length = score(childB)

            if childA_tour_length <= childB_tour_length:
                new_pop.append(childA)
            else:
                new_pop.append(childB)

        lattice.vectors = mixPops(lattice.vectors, applyMutations(new_pop, mutation_chance), basis)

        top_fitness, top_tour = testPopulation(lattice.vectors, top_fitness, top_tour)
        print(top_tour.norm())
        c_gen += 1

        tour_time_taken = time.time() - start_gen_time

    # print(tour_time_taken, time.time() - start_time + tour_time_taken)
    return top_tour, score(top_tour)


def main():
    path = "latticeBasis.txt"
    basis = Basis(readBasisFile(path))

    tour, length = runTraining(58, basis, normal_mutation_const, population_size)
    return tour, length


tour, tour_length = main()
