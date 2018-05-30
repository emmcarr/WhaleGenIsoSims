#!/usr/bin/env python

import simuPOP
import random


size = [1000,1000]
maxAge = 80
minMatingAge = 6
maxMatingAge = 60
gen = 100
nb_loci = [100]
frequency_file = "3.txt"

def runSimulation(frequency_file, sub_population_size, maxAge, minMatingAge, maxMatingAge, gen):
    '''
    sub_population_size   A vector giving the population sizes for each sub-population
    maxAge                maximum age. individuals with age > maxAge will die.
    minMatingAge          minimal mating age.
    maxMatingAge          maximal mating age.
    gen                   generations to simulate
    '''

    # Read the haplotype frequencies. There's a bit to unpack here
    # We read the lines into an array, and for each one, call split() on it to get one element per column.
    # However, we do not want this - we want the transpose, where haplotype_frequencies[0] is a vector of
    # all the frequencies for population 0, and haplotype_frequencies[1] is the corresponding vector for
    # population 2. list(map(list, zip(*t))) will achieve this transformation for us.
    # While we are at it, we also convert the strings into floats.
    with open(frequency_file, "r") as fd:
        haplotype_frequencies = list(map(list, zip(*[list(map(float, line[0:-1].split())) for line in fd])))

    if len(haplotype_frequencies) != len(sub_population_size):
        raise ValueError('The number of populations in the population size vector and the number of populations deduced from the haplotype file are different')

    sub_population_count = len(sub_population_size)
    print()
    print(sub_population_count, "subpopulations detected")

    # Now we can create the population. We want to give each population a population name, starting from A
    sub_population_names = list(map(chr, range(65, 65+sub_population_count)))
    # FIXME: Can subPopNames be a tuple here?
    pop = simuPOP.Population(sub_population_size, loci=nb_loci, infoFields=['age'], subPopNames = sub_population_names)
    sub_population_names = tuple(sub_population_names)

    # Create an attribute on each individual called 'age'. Set it to a random number between 0 and maxAge
    # Note that size is a vector - the size of each population. We have to sum these to get the total number of individuals
    individual_count = sum(sub_population_size)
    pop.setIndInfo([random.randint(0, maxAge) for x in range(individual_count)], 'age')

    # Currently we have these virtual subpopulations:
    # age < minMatingAge
    # age >= minMatingAge and age < maxMatingAge + 0.1 (age <= maxMatingAge)
    # age >= maxMatingAge + 0.1 and age < maxAge + 0.1 (maxMatingAge < age <= maxAge)
    # age >= maxAge + 0.1 (age > maxAge)
    #
    # Ideally we would want something like this:
    # 1) Immature
    # 2) Receptive female (every 3 years)
    # 3) Non-receptive female
    # 4) Mature male
    # 5) Geriatric
    # 6) Dead

    # Note that we use a cutoff InfoSplitter here, it is also possible to
    # provide a list of values, each corresponding to a virtual subpopulation.
    # FIXME: Can we call the virtual sub populations more intuitive names than 0,1,2,3?
    pop.setVirtualSplitter(simuPOP.InfoSplitter('age',
        cutoff=[minMatingAge, maxMatingAge + 0.1, maxAge + 0.1]))

    pop.evolve(
        initOps = [simuPOP.InitSex()] + [simuPOP.InitGenotype(subPops = sub_population_names[i], freq=haplotype_frequencies[i]) for i in range(0, sub_population_count)],
        # increase age by 1
        preOps = simuPOP.InfoExec('age += 1'),
        matingScheme = simuPOP.HeteroMating(
            # age <= maxAge, copy to the next generation (weight=-1)
            # subPops is a list of tuples that will participate in mating. The tuple is a pair (subPopulation, virtualSubPopulation)
            # First, we propagate (clone) all individuals in all subpopulations (and all VSPs except the ones who are now in the VSP of deceased individuals) to the next generation
            [simuPOP.CloneMating(subPops=[(sub_population, virtual_sub_population) for sub_population in range(0, sub_population_count) for virtual_sub_population in [0,1,2]], weight=-1),
            # Then we simulate random mating only in VSP 1 (ie sexually active individuals)
            simuPOP.RandomMating(subPops=[(sub_population, 1) for sub_population in range(0, sub_population_count)], weight=1)]),
        postOps = [
            # count the individuals in each virtual subpopulation
            # simuPOP.Stat(popSize=True, subPops=[(0,0), (0,1), (0,2), (0,3), (1,0), (1, 1), (1, 2), (1, 3)]),
            # print virtual subpopulation sizes (there is no individual with age > maxAge after mating)
            # simuPOP.PyEval(r"'Size of age groups: %s\n' % (','.join(['%d' % x for x in subPopSize]))")

            # Alternatively, calculate the Fst
            # FIXME: How does this actually work? Does it work for > 2 populations? I don't really understand it yet
            simuPOP.Stat(structure=range(1), subPops=sub_population_names, suffix='_AB', step=10),
            simuPOP.PyEval(r"'Fst=%.3f \n' % (F_st_AB)", step=10)

        ],
        gen = gen
    )

if __name__ == '__main__':
    runSimulation(frequency_file, size, maxAge, minMatingAge, maxMatingAge, gen)
