#!/usr/bin/env python

import simuPOP
import random
import numpy

size = [100,100]
minMatingAge = 6
maxMatingAge = 50
years = 200
nb_loci = 100
scenario_id = "1"

mean_C = [16.7, 20.5]
variance_C = [3.24, 2.89]
mean_N = [5.5, 8.7]
variance_N = [0.25, 0.49]

numberOfFeedingGrounds = 2

def computeRatios(pop):
    i = 0
    for i in range(0, pop.numSubPop()):
        for indiv in pop.individuals(i):
            # The 'feeding_ground' info field is a float. We cannot use that as an array index. Convert to an int
            feeding_ground = int(indiv.info('feeding_ground'))
            indiv.setInfo(numpy.random.normal(mean_C[feeding_ground], numpy.sqrt(variance_C[feeding_ground])), 'carbon')
            indiv.setInfo(numpy.random.normal(mean_N[feeding_ground], numpy.sqrt(variance_N[feeding_ground])), 'nitrogen')
    return True


def runSimulation(scenario_id, sub_population_size, minMatingAge, maxMatingAge, gen):
    '''
    sub_population_size   A vector giving the population sizes for each sub-population
    minMatingAge          minimal mating age.
    maxMatingAge          maximal mating age. Individuals older than this are effectively dead
    years                 number of years to simulate
    '''

    # scenario_id describes the batch of files to load
    # The mitochondrial DNA will be in mtdna_<scenario_id>
    # The SNP DNA will be in snp_<scenario_id>

    # Read the mitochondrial haplotype frequencies. There's a bit to unpack here
    # We read the lines into an array, and for each one, call split() on it to get one element per column.
    # However, we do not want this - we want the transpose, where haplotype_frequencies[0] is a vector of
    # all the frequencies for population 0, and haplotype_frequencies[1] is the corresponding vector for
    # population 2. list(map(list, zip(*t))) will achieve this transformation for us.
    # While we are at it, we also convert the strings into floats.
    mitochondrial_file = "mtdna_" + scenario_id + ".txt"
    with open(mitochondrial_file, "r") as fd:
        haplotype_frequencies = list(map(list, zip(*[list(map(float, line[0:-1].split())) for line in fd])))

    if len(haplotype_frequencies) != len(sub_population_size):
        raise ValueError('The number of populations in the population size vector and the number of populations deduced from the haplotype file are different')

    # Now read the SNP data. This builds a 2D array indexed as snp[locus][population]
    snp_file = "snp_" + scenario_id + ".txt"
    with open(snp_file, "r") as fd:
        snp = [list(map(float, line[0:-1].split())) for line in fd]

    sub_population_count = len(sub_population_size)
    print()
    print(sub_population_count, "subpopulations detected")

    # Now we can create the population. We want to give each population a population name, starting from A
    sub_population_names = list(map(chr, range(65, 65+sub_population_count)))
    # FIXME: Can subPopNames be a tuple here? ELC: could we set number of sub_pops as a global variable?
    # We have two chromosomes. The first is an autosome with nb_loci loci, and the second is the mitochondrial chromosome with 1 locus
    pop = simuPOP.Population(sub_population_size,
                                 ploidy=2,
                                 loci=[nb_loci, 1],
                                 ancGen=2,
                                 infoFields=['age', 'ind_id', 'father_id', 'mother_id', 'nitrogen', 'carbon', 'feeding_ground'],
                                 subPopNames = sub_population_names,
                                 chromTypes=[simuPOP.AUTOSOME, simuPOP.MITOCHONDRIAL])
    sub_population_names = tuple(sub_population_names)

    # Create an attribute on each individual called 'age'. Set it to a random number between 0 and maxMatingAge
    # Note that size is a vector - the size of each population. We have to sum these to get the total number of individuals
    individual_count = sum(sub_population_size)

    # Assign a random age to each individual
    pop.setIndInfo([random.randint(0, maxMatingAge) for x in range(individual_count)], 'age')
    # Assign a random feeding ground to each individual
    pop.setIndInfo([random.randint(0, numberOfFeedingGrounds-1) for x in range(individual_count)], 'feeding_ground')


    # Currently we have these virtual subpopulations:
    # age < minMatingAge (juvenile)
    # age >= minMatingAge and age < maxMatingAge + 0.1 (age <= maxMatingAge) (mature)
    # age >= maxMatingAge (dead)
    #
    # Ideally we would want something like this:
    # 1) Immature
    # 2) Receptive female (every 3 years)
    # 3) Non-receptive female
    # 4) Mature male
    # 5) Dead
    # ELC: but perhaps we can just give females a 1/3 chance of reproducing instead?
    # Note that we use a cutoff InfoSplitter here, it is also possible to
    # provide a list of values, each corresponding to a virtual subpopulation.
    # FIXME: Can we call the virtual sub populations more intuitive names than 0,1,2,3?
    pop.setVirtualSplitter(simuPOP.InfoSplitter('age',
        cutoff=[minMatingAge, maxMatingAge + 0.1]))


    pop.evolve(
        initOps = [simuPOP.InitSex(), simuPOP.IdTagger()] +
                       [simuPOP.InitGenotype(subPops = sub_population_names[i], freq=haplotype_frequencies[i], loci=[nb_loci]) for i in range(0, sub_population_count)] +
                       [simuPOP.InitGenotype(subPops = sub_population_names[i], freq=[snp[n][i], 1-snp[n][i]], loci=[n]) for i in range(0, sub_population_count) for n in range(0, nb_loci-1)],
        # increase age by 1
        preOps = simuPOP.InfoExec('age += 1'),
        matingScheme = simuPOP.HeteroMating(
            # age <= maxAge, copy to the next generation (weight=-1)
            # subPops is a list of tuples that will participate in mating. The tuple is a pair (subPopulation, virtualSubPopulation)
            # First, we propagate (clone) all individuals in all subpopulations (and all VSPs except the ones who are now in the VSP of deceased individuals) to the next generation
            [simuPOP.CloneMating(ops=[simuPOP.CloneGenoTransmitter(chroms=[0,1])],
                                     subPops=[(sub_population, virtual_sub_population) for sub_population in range(0, sub_population_count) for virtual_sub_population in [0,1]], weight=-1),
            # Then we simulate random mating only in VSP 1 (ie reproductively mature individuals)
            simuPOP.RandomMating(ops=[simuPOP.MitochondrialGenoTransmitter(),
                                      simuPOP.MendelianGenoTransmitter(),
                                      simuPOP.IdTagger(),
                                      simuPOP.InheritTagger(mode=simuPOP.MATERNAL, infoFields=['feeding_ground']),
                                      simuPOP.PedigreeTagger()],
                                 subPops=[(sub_population, 1) for sub_population in range(0, sub_population_count)], weight=1)]),
        postOps = [

        # Determine the isotopic ratios in individuals
        simuPOP.PyOperator(func=computeRatios)

            # count the individuals in each virtual subpopulation
            #simuPOP.Stat(popSize=True, subPops=[(0,0), (0,1), (0,2), (1,0), (1, 1), (1, 2)]),
            # print virtual subpopulation sizes (there is no individual with age > maxAge after mating)
            #simuPOP.PyEval(r"'Size of age groups: %s\n' % (','.join(['%d' % x for x in subPopSize]))")

            # Alternatively, calculate the Fst
            # FIXME: How does this actually work? Does it work for > 2 populations? I don't really understand it yet
            # ELC: it is a calculation that partitions variance among and between populations, and can be calculated as a 
            # global statistic or on a pairwise basis. We use it as an indication of genetic differentiation.

#            simuPOP.Dumper(structure=False),
#            simuPOP.Stat(structure=range(1), subPops=sub_population_names, suffix='_AB', step=10),
#            simuPOP.PyEval(r"'Fst=%.3f \n' % (F_st_AB)", step=10)
        ],
        gen = years
    )

    ped = simuPOP.Pedigree(pop);
    print("This is the pedigree stuff")
    simuPOP.dump(pop);

    return;



if __name__ == '__main__':
    runSimulation(scenario_id, size, minMatingAge, maxMatingAge, years)


# Plan
# * Add simulation of mitochondrial DNA
# * Break populations into VSPs and assign breeding scheme amongst VSPs
# * SNP


# Subpopulation 1: individuals who breed at wintering ground 1
# Subpopulation 2: individuals who breed at wintering ground 2

# VSP 1: individuals who feed at feeding ground 1
# VSP 2: individuals who feed at feeding ground 2



#   Migration ***
# Only males migrate
# Migrants spend 1 year in the opposing breeding grounds then return to their home population
# On average 1% of the males migrate to the other breeding ground in a given season


#   Isotopic information
# Each individual has 3 info values:
#   * Feeding ground (F1 or F2)
#   * Carbon value
#   * Nitrogen value
#  The feeding ground is fixed at birth an inherited from the mother
#  The carbon and nitrogen values are sampled at birth from the feeding ground of the mother. Each feeding ground is associated with a scalar distribution of these two values.
