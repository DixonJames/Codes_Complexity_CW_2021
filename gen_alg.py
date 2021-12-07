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



def genRandomTour(city_num):
    tour = [num for num in range(int(city_num))]
    random.shuffle(tour)
    return tour

def genStartPopulation(population_num, basis:Basis):
    population = []
    for p in range(population_num):
        rand_tour = basis.genRandPoint()
        population.append(rand_tour)
    return population

def tourFitness(v:LatticeVector):
    return v.norm()

def findDuplicate(list):
    dupes = []
    for char in list:
        if list.count(char) != 1:
            dupes.append(char)
    return dupes

def findMissing(original, trial):
    missing = []
    for char in original:
        if char not in trial:
            missing.append(char)
    return list(set(missing))


def basicCrossoverTours(lattice, p, q):
    return lattice.modSivAvg(p, q), lattice.sivDiff(p, q)


def basicSwapMutation(v:LatticeVector):
    A_basis_vector = random.randint(0, len(v.coeffs)-1)
    B_basis_vector = random.randint(0, len(v.coeffs)-1)

    while B_basis_vector == A_basis_vector:
        B_basis_vector = random.randint(0, len(v.coeffs) - 1)


    v.coeffs[A_basis_vector], v.coeffs[B_basis_vector] = v.coeffs[B_basis_vector], v.coeffs[A_basis_vector]

    return v



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
        if top_fitness >= tourFitness(tour):
            top_fitness = tourFitness(tour)
            top_tour = tour
    return top_fitness, top_tour


def runTraining(time_frame, basis:Basis,  mutation_chance, popsize):
    start_time = time.time()
    tau = LatticeVector(basis, [9999 for i in range(basis.dim)]).norm()
    top_fitness = LatticeVector(basis, [9999 for i in range(basis.dim)]).norm()
    top_tour = []
    population = genStartPopulation(popsize, basis)

    lattice = Lattice(basis)
    lattice.vectors = population


    elapsed_time = time.time() - start_time
    c_gen = 0

    tour_time_taken = time.time() - time.time()

    while (time.time() - start_time + tour_time_taken < time_frame):
        start_gen_time = time.time()


        population_fitness = [tourFitness(v) for v in lattice.vectors]
        total_fitness = sum(population_fitness)
        population_percentage = [(tau - fitness)/total_fitness for fitness in population_fitness]

        new_pop = []

        for j in range(popsize):
            #making as many children as there are parents
            parents = random.choices(lattice.vectors, weights=population_percentage, k=2)
            childA, childB = basicCrossoverTours(lattice, parents[0], parents[1])

            childA_tour_length = tourFitness(childA)
            childB_tour_length = tourFitness(childB)

            if childA_tour_length <= childB_tour_length:
                new_pop.append(childA)
            else:
                new_pop.append(childB)

        lattice.vectors = applyMutations(new_pop, mutation_chance)


        top_fitness, top_tour = testPopulation(lattice.vectors, top_fitness, top_tour)
        print(top_tour.norm())
        c_gen += 1

        tour_time_taken = time.time() - start_gen_time

    #print(tour_time_taken, time.time() - start_time + tour_time_taken)
    return top_tour, tourFitness(top_tour)

def main():
    path = "latticeBasis.txt"
    basis = Basis(readBasisFile(path))

    tour, length = runTraining(58, basis, normal_mutation_const, population_size)
    return tour, length



tour, tour_length = main()