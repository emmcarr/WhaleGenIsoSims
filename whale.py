#!/usr/bin/env python

## this code is to generate simulated whale populations that have mtDNA, SNP genotypes and stable isotope profiles
## the code generates two subpopulations in simupop that correspond to whale wintering grounds
## virtual subpopulations (VSPs) are used to represent feeding ground and age-sex classes

import simuPOP
from simuPOP.sampling import drawRandomSample
from simuPOP.utils import export
import random
import numpy
from operator import add

## set parameters for simulation - breeding ground which are subpopulations in simupop terms
size = [100,100] ## size of each of the two subpopulations
minMatingAge = 6 ## minimum age at first reproduction
maxMatingAge = 50 ## max age of reproduction
years = 60 
nb_loci = 100 ## number of loci to simulate
scenario_id = "1"

## setting up the feeding ground variables
## mean (mean_) and variance (variance_) set for both C and N for two feeding grounds
## deviant proportion: proportion of males that will go to non-natal wintering ground for one winter/breeding opportunity
mean_C = [16.7, 20.5] 
variance_C = [3.24, 2.89]
mean_N = [5.5, 8.7]
variance_N = [0.25, 0.49]
deviant_proportion = 0.1

## Sample count is number of samples taken per wintering ground
numberOfFeedingGrounds = 2
sample_count = 60

# Needed a way to ensure that the simulations begin with whales that are related and show correlation between 
# SI and genetic data. Did this by rapidly expanding population, creating many offspring
# For the first 10 generations, we expand the next generation by 7% (this leads to a rough doubling after 10 years).
# After the population_growing_period, each subpopulation size is kept constant
population_growth_rate = 1.07
population_growing_period = 10 # in years

def postop_processing(pop):
    for i in range(0, pop.numSubPop()):
        for individual in pop.individuals(i):
            # The feeding ground is fixed at birth (inherited from mother)
            # The C and N values are sampled from a distribution based on the feeding ground each year
            # The 'feeding_ground' info field is a float. We cannot use that as an array index so convert to an int
            if (years == 1):
                feeding_ground = int(individual.info('feeding_ground'))
                individual.setInfo(numpy.random.normal(mean_C[feeding_ground], numpy.sqrt(variance_C[feeding_ground])), 'carbon')
                individual.setInfo(numpy.random.normal(mean_N[feeding_ground], numpy.sqrt(variance_N[feeding_ground])), 'nitrogen')
            else:
                simuPOP.InheritTagger(mode=simuPOP.MATERNAL, infoFields=['nitrogen', 'carbon', 'native_breeding_ground', 'feeding_ground'])
            # print("Individual ", individual.info('ind_id'), " has native breeding ground ", individual.info('native_breeding_ground'), " and is currently at breeding ground ", i)
            # Migration
            # Initially, set the migrate_to to the current population of the individual
            individual.setInfo(i, 'migrate_to')
            # If the individual is a male, then we can optionally migrate them using the migrate_to info field
            if individual.sex() == simuPOP.MALE and individual.info('age') >= minMatingAge:
                # If the individual has already migrated, always move them back
                if individual.info('native_breeding_ground') != i:
                    #print("Moving individual ", individual.info('ind_id'), " back to their native breeding ground ", individual.info('native_breeding_ground'), " from temporary breeding ground ", i)
                    individual.setInfo(individual.info('native_breeding_ground'), 'migrate_to')
                # Otherwise, migrate them to another population with a probabilistic model
                elif numpy.random.uniform() < deviant_proportion:
                    # Individual will migrate.
                    new_population = (i + 1) % 2
                    #print("Individual ", individual.info('ind_id'), " will migrate to ", new_population)
                    individual.setInfo(new_population, 'migrate_to')
    return True

def init_native_breeding_grounds(pop):
    # Assign the native breeding ground to each individual. I don't know how to do this except by doing it individually
    # Fortunately, we can just inherit this maternally, so it only has to be run once
    for i in range(0, pop.numSubPop()):
        for individual in pop.individuals(i):
            individual.setInfo(i, 'native_breeding_ground');
    return True

def configure_new_population_size(gen, pop):
    # It is critical to specify the sub population sizes independently of each other. Each sub-population may be a different size
    if (gen < population_growing_period):
        return [pop.subPopSize(0) * population_growth_rate, pop.subPopSize(1) * population_growth_rate]
    else:
        return [pop.subPopSize(0), pop.subPopSize(1)]

def runSimulation(scenario_id, sub_population_size, minMatingAge, maxMatingAge, gen):
    '''
    sub_population_size   A vector giving the population sizes for each sub-population. The subpopulations determine which breeding ground an individual belongs to
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
    
    # We have two chromosomes. The first is an autosome with nb_loci loci, and the second is the mitochondrial chromosome with 1 locus
    pop = simuPOP.Population(sub_population_size,
                                 ploidy=2,
                                 loci=[nb_loci, 1],
                                 ancGen=2,
                                 infoFields=['age', 'ind_id', 'father_id', 'mother_id', 'nitrogen', 'carbon', 'feeding_ground', 'native_breeding_ground', 'migrate_to'],
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
    # 
    # Note that we use a cutoff InfoSplitter here, it is also possible to
    # provide a list of values, each corresponding to a virtual subpopulation.
    pop.setVirtualSplitter(simuPOP.CombinedSplitter([
        simuPOP.ProductSplitter([simuPOP.SexSplitter(),
                                 simuPOP.InfoSplitter('age', cutoff=[minMatingAge, maxMatingAge + 0.1], names=['juvenile', 'mature', 'dead'])])],
                                 vspMap = [[0], [1], [2], [3], [4], [5], [0, 1, 3, 4], [1,4]],
                                  names = ['Juvenile Male', 'Mature Male', 'Dead Male', 'Juvenile Female', 'Mature Female', 'Dead Female', 'Not dead yet', 'Active']))


    pop.evolve(
        initOps = [simuPOP.InitSex(), simuPOP.IdTagger(), simuPOP.PyOperator(func=init_native_breeding_grounds)] +
                       [simuPOP.InitGenotype(subPops = sub_population_names[i], freq=haplotype_frequencies[i], loci=[nb_loci]) for i in range(0, sub_population_count)] +
                       [simuPOP.InitGenotype(subPops = sub_population_names[i], freq=[snp[n][i], 1-snp[n][i]], loci=[n]) for i in range(0, sub_population_count) for n in range(0, nb_loci-1)],
        # increase age by 1
        preOps = [simuPOP.InfoExec('age += 1')],
        matingScheme = simuPOP.HeteroMating([
            # age <= maxAge, copy to the next generation (weight=-1)
            # subPops is a list of tuples that will participate in mating. The tuple is a pair (subPopulation, virtualSubPopulation)
            # First, we propagate (clone) all individuals in all subpopulations (and all VSPs except the ones who are now in the VSP of deceased individuals) to the next generation
            simuPOP.CloneMating(ops=[simuPOP.CloneGenoTransmitter(chroms=[0,1])],
                                subPops=[(sub_population, 6) for sub_population in range(0, sub_population_count)],
                                weight=-1),
            # Then we simulate random mating only in VSP 1 (ie reproductively mature individuals) within subpopulation (breeding/winter grounds)
            simuPOP.RandomMating(ops=[simuPOP.MitochondrialGenoTransmitter(),
                                      simuPOP.MendelianGenoTransmitter(),
                                      simuPOP.IdTagger(),
                                      simuPOP.InheritTagger(mode=simuPOP.MATERNAL, infoFields=['feeding_ground']),
                                      simuPOP.InheritTagger(mode=simuPOP.MATERNAL, infoFields=['native_breeding_ground']),
                                      simuPOP.PedigreeTagger()],
                                 subPops=[(sub_population, 7) for sub_population in range(0, sub_population_count)],
                                 weight=1)],
                subPopSize=configure_new_population_size
                                                ),
        postOps = [

        # Determine the isotopic ratios in individuals
        simuPOP.PyOperator(func=postop_processing),
        simuPOP.Migrator(mode=simuPOP.BY_IND_INFO),
            # count the individuals in each virtual subpopulation
            #simuPOP.Stat(popSize=True, subPops=[(0,0), (0,1), (0,2), (1,0), (1, 1), (1, 2)]),
            # print virtual subpopulation sizes (there is no individual with age > maxAge after mating)
            #simuPOP.PyEval(r"'Size of age groups: %s\n' % (','.join(['%d' % x for x in subPopSize]))")

            # Alternatively, calculate the Fst
            # FIXME: How does this actually work? Does it work for > 2 populations? I don't really understand it yet
            # ELC: it is a calculation that partitions variance among and between populations, and can be calculated as a 
            # global statistic or on a pairwise basis. We use it as an indication of genetic differentiation.

            simuPOP.Stat(structure=range(1), subPops=sub_population_names, suffix='_AB', step=10),
            simuPOP.PyEval(r"'Fst=%.3f \n' % (F_st_AB)", step=10)
        ],
        gen = years
    )

    #simuPOP.dump(pop, width=3, loci=[], subPops=[(simuPOP.ALL_AVAIL, simuPOP.ALL_AVAIL)], max=1000, structure=False);
    #return



    ped = simuPOP.Pedigree(pop);
    print("This is the pedigree stuff")
    simuPOP.dump(pop);

    # Now sample the individuals
    sample = drawRandomSample(pop, sizes=[sample_count]*sub_population_count)

    # Print out the allele frequency data
    simuPOP.stat(sample, alleleFreq=simuPOP.ALL_AVAIL)
    frequencies = sample.dvars().alleleFreq;
    with open('freq.txt', 'w') as freqfile:
        index = 0
        for locus in frequencies:
            if (locus == nb_loci):
                continue
            if (len(frequencies[locus]) < 2):
                continue
            print(index, end=' ', file=freqfile)
            index = index + 1
            for allele in frequencies[locus]:
                print(frequencies[locus][allele], end=' ', file=freqfile)
            print(file=freqfile)

    # We want to remove monoallelic loci. This means a position in the genotype for which all individuals have the same value in both alleles
    # To implement this we will build up a list of loci that get ignored when we dump out the file. Generally speaking, if we add all the values up
    # then either they will sum to 0 (if all individuals have type 0) or to the number of individuals * 2 (if all individuals have type 1)
    geno_sum = [0] * (nb_loci + 1) * 2;
    for individual in sample.individuals():
        geno_sum = list(map(add, geno_sum, individual.genotype()))
    final_sum = list(map(add, geno_sum[:(nb_loci+1)], geno_sum[(nb_loci+1):]))

    monoallelic_loci = [];
    for i in range(0, nb_loci):
        if final_sum[i] == 0 or final_sum[i] == sample_count*sub_population_count*2:
            monoallelic_loci = [i] + monoallelic_loci
    monoallelic_loci = sorted(monoallelic_loci, reverse=True)

    nb_ignored_loci = len(monoallelic_loci)
    # Generate the two files
    with open('mixfile.txt', 'w') as mixfile:
        with open('haploiso.txt', 'w') as haplofile:
            print(sub_population_count, nb_loci - nb_ignored_loci, 2, 1, file=mixfile)
            print("sex, haplotype, iso1, iso2, native_ground", file=haplofile);
            for i in range(0, nb_loci - nb_ignored_loci):
                print('Loc', i+1, sep='_', file=mixfile);
            for individual in sample.individuals():
                genotype = individual.genotype();
                print(1 if individual.sex() == 1 else 0,
                      genotype[nb_loci],
                      individual.info('carbon'),
                      individual.info('nitrogen'),
#                      int(individual.info('native_breeding_ground')),
                      file=haplofile, sep=' ')
                print(int(individual.info('native_breeding_ground')+1), end=' ', file=mixfile)
                for i in range(0, nb_loci):
                    if i not in monoallelic_loci:
                        print(genotype[i]+1, genotype[i+nb_loci+1]+1, ' ', end='', sep='', file=mixfile)
                print(file=mixfile);
    return sample




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

