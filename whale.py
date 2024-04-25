#This script uses simuPOP to simulate the evolution of Aotearoa New Zealand's
#population of tohorā Southern Right Whales and its recovery from a very
#harsh bottleneck.
#Code originally written by Dr Emma Carroll (EC) and later edited by
#Dr Marina Klemm (MK). This code also uses some chunks shared by other
#simuPOP users, cited where needed.

#Contact: marinaklemm@gmail.com
#Github: marina-klemm


# =============================================================================
# Loading all packages
# =============================================================================
import os
import csv
import string
import random
import logging
from operator import add
from pathlib import Path

import numpy as np
import matplotlib

import simuPOP as sim
from simuPOP.sampling import drawRandomSample
from simuPOP.utils import export
from simuPOP.utils import Trajectory, simulateForwardTrajectory, export, Exporter
from simuPOP.utils import viewVars

import simuOpt
simuOpt.setOptions(quiet=True)

# get the random seed from simuPOP
# create filenames that don't collide (i.e. overwrite each other)
SEED = "0x%08x" % sim.getRNG().seed()

# WHERE to store results
RESULTS_DIR = Path("./results")
if not RESULTS_DIR.exists():
    RESULTS_DIR.mkdir()



# =============================================================================
# MK: Population numbers used as a reference here
# =============================================================================
# Used Jackson et al, 2016 to estimate trajectory in 25 years intervals
# Starting in 1830 and ending in 2030: (divided into 25 years)
    
#1830-1855: ~ 30,000 individuals
#1855-1880: ~ 3,000
#1880-1905: ~ 200
#1905-1930: ~ 150
#1930-1955: ~ 40
#1955-1980: ~ 150
#1980-2005: ~ 1000 
#2005-2030: ~ 2,139 (2009)
#2010-2030: ~ 2500 (maybe, at least?)

   
# I will start to print them at generation 75, so the first 74 generations
# will be a "burn in" time.

#Gen 75: 1830-1855: ~ 30,000 individuals
#Gen 76: 1855-1880: ~ 3,000
#Gen 77: 1880-1905: ~ 200
#Gen 78: 1905-1930: ~ 150
#Gen 79: 1930-1955: ~ 40
#Gen 80: 1955-1980: ~ 150
#Gen 81: 1980-2005: ~ 1000 
#Gen 82: 2005-2030: ~ 2,139 (2009)
#Gen 83: 2030-2055: ~ 2500 (maybe, hopefully, at least?)


   
#Using demo function from simuBottlenecks_v2_serial.py by Juhana Kamonnen

staticPhaseEnd = 75 #year 1830, before whaling
firstDecline	= 76 #year 1855, whaling starts
secondDecline	= 77 #year 1880, whaling continues
thirdDecline	= 78 #year 1905, whaling continues
bottleNeck = 79 #year 1930, tohorā almost extinct
firstRecovery	= 80 #year 1955, whales start recovering
secondRecovery	= 81 #year 1980, there is hope
currentGen	= 82 #year 2005, current population doubles from last generation
hopefulFuture	= 83 #year 2030, can we be this hopeful?


#Population sizes
pop_staticPhaseEnd = 30000 #year 1830, before whaling
pop_firstDecline	= 3000 #year 1855, whaling starts
pop_secondDecline	= 200 #year 1880, whaling continues
pop_thirdDecline	= 150 #year 1905, whaling continues
pop_bottleNeck = 40  #year 1930, tohorā almost extinct
pop_firstRecovery	= 150 #year 1955, whales start recovering
pop_secondRecovery	= 1000 #year 1980, there is hope
pop_currentGen	= 2138 #year 2005, current population doubles from last generation
pop_hopefulFuture	= 4000 #year 2030, can we be this hopeful?


def get_filename(filename, seed=SEED, results_dir=RESULTS_DIR):
    """
    Return a unique filename in the output directory with the random seed preprepended:
    
    >>> get_filename("test.csv")
    results/a26d85ecf7d490bc9e44eeb8c80f00d8_test.csv

    >>> get_filename("result.dat")
    results/a26d85ecf7d490bc9e44eeb8c80f00d8_result.dat
    """
    return Path(results_dir) / ("%s_%s" % (seed, filename))


def demo(gen, pop):  # demographics function to control population growth (actually birth rate now).
                     # number of bottlenecks and the severities are as agreed on 01.09.2009. 
                     #Check the generation "names" from the beginning of the script

    if gen < staticPhaseEnd:
        return [pop_staticPhaseEnd, pop_staticPhaseEnd]
    
    elif gen < firstDecline:
        return [pop_firstDecline, pop_firstDecline] 
    
    elif gen < secondDecline:
        return [pop_secondDecline, pop_secondDecline]
    
    elif gen < thirdDecline:
        return [pop_thirdDecline, pop_thirdDecline]

    elif gen < bottleNeck:
        return [pop_bottleNeck, pop_bottleNeck]
    
    elif gen < firstRecovery:
        return [pop_firstRecovery, pop_firstRecovery]

    elif gen < secondRecovery:
        return [pop_secondRecovery, pop_secondRecovery]

    elif gen < currentGen:
        return [pop_currentGen, pop_currentGen]
    
    else:
        return [pop_hopefulFuture, pop_hopefulFuture]


# =============================================================================
# I need to add Juhana Kammonen's overspill function too, otherwise I cannot add
# weight=-1 to CloneMating because the population cannot handle that number of
# individuals (since there are some generations with a significant reduction in
# the number of individuals)
# Error message:   
# ValueError: mating.cpp: 1886 Mating scheme with a negative weight of 1 would 
# like to produce 14696 offspring, but there are only 1500 unclaimed offspring left.
# This error has now been fixed with Bo Peng's help and a new patch #issue114.
# =============================================================================

#Ok, after adding the lethalEvent, overspillLethalEvent and removeOverspill,
#I stopped having the 1886 error (above) but started finding the 1906 error:
#1906 An exact (all negative) weight system is used, but does not fill offspring subpopulation.

#Following what was recommended by Bo in the simuPOP issue #75,
#https://github.com/BoPeng/simuPOP/issues/75,
#I will add weight=0 to the second mating scheme.


def overspillLethalEvent(pop, numRemovables):

    if numRemovables <= 0:
        return True
    else:
        numRemovables = numRemovables + 1
        print('Removables: %d' % numRemovables)
        indices = range(pop.popSize()) # a list of individual indices
        print('Indices: %d' % len(indices))
        affected = random.sample(indices, int(numRemovables))
        pop.removeIndividuals(affected)
    return True



def removeOverspill(pop):

    gen = pop.dvars().gen
    
    if gen <= staticPhaseEnd:
        nextSubPopSizes = [pop_staticPhaseEnd, pop_staticPhaseEnd]
        maxAllowed = sum(nextSubPopSizes)
        numRemovables = pop.popSize() - maxAllowed 
        overspillLethalEvent(pop, numRemovables)

    elif gen == firstDecline:
        nextSubPopSizes = [pop_firstDecline, pop_firstDecline]  
        maxAllowed = sum(nextSubPopSizes)
        numRemovables = pop.popSize() - maxAllowed 
        overspillLethalEvent(pop, numRemovables)
    
    elif gen == secondDecline:
        nextSubPopSizes = [pop_secondDecline, pop_secondDecline]
        maxAllowed = sum(nextSubPopSizes)
        numRemovables = pop.popSize() - maxAllowed
        overspillLethalEvent(pop, numRemovables)
        
    elif gen == thirdDecline:
        nextSubPopSizes = [pop_thirdDecline, pop_thirdDecline]
        maxAllowed = sum(nextSubPopSizes)
        numRemovables = pop.popSize() - maxAllowed
        overspillLethalEvent(pop, numRemovables)
        
    elif gen == bottleNeck:
        nextSubPopSizes = [pop_bottleNeck, pop_bottleNeck]
        maxAllowed = sum(nextSubPopSizes)
        numRemovables = pop.popSize() - maxAllowed
        overspillLethalEvent(pop, numRemovables)
        
    elif gen == firstRecovery:
        nextSubPopSizes = [pop_firstRecovery, pop_firstRecovery]
        maxAllowed = sum(nextSubPopSizes)
        numRemovables = pop.popSize() - maxAllowed
        overspillLethalEvent(pop, numRemovables)
        
    elif gen == secondRecovery:
        nextSubPopSizes = [pop_secondRecovery, pop_secondRecovery]
        maxAllowed = sum(nextSubPopSizes)
        numRemovables = pop.popSize() - maxAllowed
        overspillLethalEvent(pop, numRemovables)
        
    elif gen == currentGen:
        nextSubPopSizes =  [pop_currentGen, pop_currentGen]
        maxAllowed = sum(nextSubPopSizes)
        numRemovables = pop.popSize() - maxAllowed
        overspillLethalEvent(pop, numRemovables)
        
    else: #final rise
        nextSubPopSizes = [pop_hopefulFuture, pop_hopefulFuture]
        maxAllowed = sum(nextSubPopSizes)
        numRemovables = pop.popSize() - maxAllowed
        overspillLethalEvent(pop, numRemovables)
        
    gen += 1
    return True

# =============================================================================
# Report mom and dad (from BO's github issue #25)
# =============================================================================
#This is just to debug it; it is added as a sim.PyOperator((report)) inside the
#preOps in the pop.evolve function, but I commented it out right now to avoid all
#the printing it does. 

def report(dad, mom, off):
   # print('{} ({}) from {} {}'.format(off.ind_id, off.sex(), dad.ind_id, mom.ind_id))
   #instead of printing it into the terminal, create a file for debugging:
   with open(get_filename('parents_check.txt'), 'a') as file:
       file.write('{} ({}) from {} {}\n'.format(off.ind_id, off.sex(), dad.ind_id, mom.ind_id))
   return True





# =============================================================================
# Initial code (by EC, edited by MK, unnused things removed)
# =============================================================================
## EC
## this code is to generate simulated whale populations that have mtDNA, SNP genotypes
## and stable isotope profiles. the code generates two subpopulations in sim that 
## correspond to whale wintering grounds 
## virtual subpopulations (VSPs) are used to represent feeding ground and age-sex classes

## set parameters for simulation - breeding ground which are subpopulations in simupop terms

minMatingAge = 6 ## minimum age at first reproduction
maxMatingAge = 50 ## max age of reproduction
gen = 82 ## for trajectory, based on the model6
nb_loci = 100 ## number of loci to simulate
scenario_id = "1"

## setting up the feeding ground variables
## mean (mean_) and variance (variance_) set for both C and N for two feeding grounds
## deviant proportion: proportion of males that will go to non-natal wintering 
# ground for one winter/breeding opportunity
mean_C = [16.7, 20.5] 
variance_C = [3.24, 2.89]
mean_N = [5.5, 8.7]
variance_N = [0.25, 0.49]
deviant_proportion = 0.1

## Sample count is number of samples taken per wintering ground
numberOfFeedingGrounds = 1
sample_count = 60


# Needed a way to ensure that the simulations begin with whales that are related 
# and show correlation between 
# SI and genetic data. Did this by rapidly expanding population, creating many offspring
# For the first 10 generations, we expand the next generation by 7% (this leads 
# to a rough doubling after 10 years).
# After the population_growing_period, each subpopulation size is kept constant
population_growth_rate = 1.07 
population_growing_period = 10 # in years


def postop_processing(pop):
    for i in range(0, pop.numSubPop()):
        for individual in pop.individuals(i):
            # The feeding ground is fixed at birth (inherited from mother)
            # The C and N values are sampled from a distribution based on the 
            # feeding ground each year
            # The 'feeding_ground' info field is a float. We cannot use that as 
            # an array index so convert to an int
            feeding_ground = int(individual.info('feeding_ground'))
            individual.setInfo(np.random.normal(mean_C[feeding_ground], np.sqrt(variance_C[feeding_ground])), 'carbon')
            individual.setInfo(np.random.normal(mean_N[feeding_ground], np.sqrt(variance_N[feeding_ground])), 'nitrogen')

            # print("Individual ", individual.info('ind_id'), " has native breeding ground ", 
            # individual.info('native_breeding_ground'), " and is currently at breeding ground ", i)
            # Migration
            # Initially, set the migrate_to to the current population of the individual
            individual.setInfo(i, 'migrate_to')
            # If the individual is a male, then we can optionally migrate them using the migrate_to info field
            if individual.sex() == sim.MALE and individual.info('age') >= minMatingAge:
                # If the individual has already migrated, always move them back
                if individual.info('native_breeding_ground') != i:
                    #print("Moving individual ", individual.info('ind_id'), " back to their native breeding ground ", individual.info('native_breeding_ground'), " from temporary breeding ground ", i)
                    individual.setInfo(individual.info('native_breeding_ground'), 'migrate_to')
                # Otherwise, migrate them to another population with a probabilistic model
                elif np.random.uniform() < deviant_proportion:
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
        return [pop.subPopSize(0) * population_growth_rate]
    else:
        return [pop.subPopSize(0)]



def runSimulation(sub_population_size, minMatingAge, maxMatingAge, gen, mitochondrial_file, snp_file):
    '''
    sub_population_size   A list giving the population sizes for each sub-population. The subpopulations
                          determine which breeding ground an individual belongs to
    minMatingAge          minimal mating age.
    maxMatingAge          maximal mating age. Individuals older than this are effectively dead
    mitochondrial_file    Path to the file containing the mitochondrial data.
                          Should be a tab delimited text file with one column per subpopulation.
    snp_file              Path to the file containing the SNP data.
                          Should be a tab delimited text file with one column per subpopulation.
    '''

    # Read the mitochondrial haplotype frequencies. There's a bit to unpack here
    # We read the lines into an array, and for each one, call split() on it to get one element per column.
    # However, we do not want this - we want the transpose, where haplotype_frequencies[0] is a vector of
    # all the frequencies for population 0, and haplotype_frequencies[1] is the corresponding vector for
    # population 2. list(map(list, zip(*t))) will achieve this transformation for us.
    # While we are at it, we also convert the strings into floats.
    with open(mitochondrial_file, "r") as fd:
        haplotype_frequencies = list(map(list, zip(*[list(map(float, line[0:-1].split())) for line in fd])))
    
    if len(haplotype_frequencies) != len(sub_population_size):
        raise ValueError(
            'The number of populations in the population size vector and the number of populations'
            ' deduced from the haplotype file are different'
        )

    # Now read the SNP data. This builds a 2D array indexed as snp[locus][population]
    with open(snp_file, "r") as fd:
        snp = [list(map(float, line[0:-1].split())) for line in fd]

    # how many sub_populations do we have?
    sub_population_count = len(sub_population_size)

    # Now we can create the population. 
    # We want to give each population a population name, starting from A
    sub_population_names = tuple(
        string.ascii_uppercase[i] for (i, s) in enumerate(sub_population_size)
    )
        
    # =============================================================================
    # Population
    # =============================================================================

    pop = sim.Population(
        size=pop_staticPhaseEnd,
        ploidy=2,
        loci=[nb_loci, 1],
        ancGen=2,#Number of the most recent ancestral generations 
        #to keep during evolution, i.e., ancGen=2 keep parental and
        #grandparental generations coexisting with the newest one.
        infoFields=['age', 'ind_id', 'father_id', 'mother_id', 'nitrogen', 'carbon', 'feeding_ground', 'native_breeding_ground', 'migrate_to'],
        subPopNames=sub_population_names,
        chromTypes=[sim.AUTOSOME, sim.MITOCHONDRIAL])

    # Create an attribute on each individual called 'age'. Set it to a random number 
    # between 0 and maxMatingAge
    # Note that size is a vector - the size of each population. We have to sum these
    #  to get the total number of individuals
    #i ndividual_count = sum(sub_population_size)
    individual_count = pop_staticPhaseEnd * 2

    # Assign a random age to each individual
    pop.setIndInfo(
        [random.randint(0, maxMatingAge) for x in range(individual_count)], 'age')
    # Assign a random feeding ground to each individual
    pop.setIndInfo(
        [random.randint(0, numberOfFeedingGrounds-1) for x in range(individual_count)],
        'feeding_ground')

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
    pop.setVirtualSplitter(
        sim.CombinedSplitter([
            sim.ProductSplitter([
                sim.SexSplitter(),
                sim.InfoSplitter(
                    "age",
                    cutoff=[minMatingAge, maxMatingAge + 0.1],
                    names=["juvenile", "mature", "dead"],
                ),
            ])
        ],
        vspMap=[[0], [1], [2], [3], [4], [5], [0, 1, 3, 4], [1, 4]],
        names=[
            "Juvenile Male",  # 0
            "Mature Male",    # 1
            "Dead Male",      # 2
            "Juvenile Female",# 3
            "Mature Female",  # 4
            "Dead Female",    # 5
            "Not dead yet",   # 0, 1, 3, 4
            "Active",         # 1, 4 = sexually active
        ])
    )

    ## sim.dump(pop)
    # =============================================================================
    # MK: Check if VSPs are initiallized properly
    # =============================================================================
    # pop.subPopName([0,0])
    # pop.subPopName([0,1])
    # pop.subPopName([1,7])
    # pop.subPopName([1,4])

    ##############################################################################
    # =============================================================================
    # Printing 10 alleles:
    # =============================================================================
    #initOps: applied before evolution
    #preOps: applied to parental population at the beginning of each life cycle
    #postOps: applied to offspring generation at the end of each life cycle
    #finalOps: aplied at the end of evolution

    #To fix the weight=-1 problem, I'll use BoPeng's solution
    pop.dvars().demo = demo
    sim.stat(pop, popSize=True)
    #BO: is important here because the expression 
    #popSize > demo(gen) needs demo in the population's namespace.


    ##January 14th, 2024
    #I successfully installed the patch now, and then I can try to make it work as
    #it was when I first asked for Bo's help:

    print("---")
    print("SIMON::")
    
    print(sub_population_count)
    print(sub_population_names)
    print("---")


    pop.evolve(
        initOps = [
            sim.InitSex(),
            sim.IdTagger(),
            sim.PyOperator(func=init_native_breeding_grounds)
        ] +
        [ sim.InitGenotype(subPops = sub_population_names[i], freq=haplotype_frequencies[i], loci=[nb_loci]) for i in range(0, sub_population_count)] +
        [ sim.InitGenotype(subPops = sub_population_names[i], freq=[snp[n][i], 1-snp[n][i]], loci=[n]) for i in range(0, sub_population_count) for n in range(0, nb_loci-1)],
        # increase age by 1
        preOps = [
            sim.InfoExec('age += 1'),
            sim.Stat(popSize=True), #print pop size in each generation
            sim.PyEval(r'"%s\n" % subPopSize'),
            # randomly reduce population size so that parent generation fits
            sim.PyOperator(func=removeOverspill),
            # Export population in each generation
            Exporter(
                format='csv',
                infoFields=('age', 'ind_id', 'father_id', 'mother_id', 'nitrogen', 'carbon', 'feeding_ground', 'native_breeding_ground', 'migrate_to'), 
                #output="!'dump_gen_%d.csv' % gen", step=1, begin=75
                output="!'%s_%%d.csv' %% gen" % get_filename('dump_gen'),
                #output="!'dump_gen_%d.csv' % gen",
                step=1, begin=75
            )
        ],
        matingScheme = sim.HeteroMating([
            # age <= maxAge, copy to the next generation (weight=-1)
            # subPops is a list of tuples that will participate in mating. The tuple is a pair (subPopulation, virtualSubPopulation)
            # First, we propagate (clone) all individuals in all subpopulations (and all VSPs except the ones who are now in the VSP of deceased individuals) to the next generation
            sim.CloneMating(
                ops=[sim.CloneGenoTransmitter(chroms=[0,1])],
                #MK:6 here is because two of the eight virtual subpopulations are deceased.
                subPops=[(sub_population, 6) for sub_population in range(0, sub_population_count)], 
                weight=-1
            ), 
                #MK: if weights are negative, they are multiplied to their parental subpopulation; 
                # For example: 
                #  if parental pop = (500, 1000), and weight = -2, 
                #  next generation pop= (1000, 2000). 
                # For weight -1, it keeps the number of individuals from the parental generation.
                #
                # ALSO: if there is a mix of negative and positive weights, the negative will be
                # processed first.
            # Then we simulate random mating only in VSP 1 (i.e. reproductively mature individuals)
            # within subpopulation (breeding/winter grounds)
            sim.RandomMating(
                ops=[sim.MitochondrialGenoTransmitter(),
                     sim.MendelianGenoTransmitter(),
                     sim.IdTagger(),
                     sim.InheritTagger(mode=sim.MATERNAL, infoFields=['feeding_ground']),
                     sim.InheritTagger(mode=sim.MATERNAL, infoFields=['native_breeding_ground']),
                     sim.PedigreeTagger(),
                     #sim.PyOperator((report))
                ],
                subPops=[(sub_population, 7) for sub_population in range(0, sub_population_count)],
                weight=0 #recommended by Bo:
                # The second mating scheme should have weight 0, and generate the rest of the offspring. 
                # If you get an error, it probably means the parental virtual subpopulation is empty.
             )],
            subPopSize=demo),
            # REMOVING TRAJECTORY TO SEE IF IT WORKS NOW:        
            # MK: we decided to keep the same weight as the mitochondrial transmitter.
            # sim.ControlledRandomMating(subPopSize=model10, freqFunc=traj.func(), weight=1)]
       
        postOps = [

        # Determine the isotopic ratios in individuals
        sim.PyOperator(func=postop_processing),
        sim.Migrator(mode=sim.BY_IND_INFO),
            # count the individuals in each virtual subpopulation
            #sim.Stat(popSize=True, subPops=[(0,0), (0,1), (0,2), (1,0), (1, 1), (1, 2)]),
            # print virtual subpopulation sizes (there is no individual with age > maxAge after mating)
            #sim.PyEval(r"'Size of age groups: %s\n' % (','.join(['%d' % x for x in subPopSize]))")
            # Alternatively, calculate the Fst
            # FIXME: How does this actually work? Does it work for > 2 populations? I don't really understand it yet
            # ELC: it is a calculation that partitions variance among and between populations, and can be calculated as a 
            # global statistic or on a pairwise basis. We use it as an indication of genetic differentiation.

            sim.Stat(structure=range(1), subPops=sub_population_names, suffix='_AB', step=10),
            sim.Stat(numOfMales=True, begin = 73, step = 1),
            sim.PyEval(r"'Fst=%.3f \n' % (F_st_AB)", step=10), #Print Fst every 10 steps
            sim.Stat(alleleFreq=[1, 2, 3, 4, 5, 6, 7, 8, 9, 100], vars=['alleleFreq_sp'], step=10), #added this now, to
            #calculate the allele frequencies in selected loci
           sim.PyEval(
              r"'%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n' % ("
               "subPop[0]['alleleFreq'][1][1], subPop[0]['alleleFreq'][2][1], subPop[0]['alleleFreq'][3][1],"
               "subPop[0]['alleleFreq'][4][1], subPop[0]['alleleFreq'][5][1], subPop[0]['alleleFreq'][6][1],"
               "subPop[0]['alleleFreq'][7][1], subPop[0]['alleleFreq'][8][1], subPop[0]['alleleFreq'][9][1],"
               "subPop[0]['alleleFreq'][100][1], subPop[1]['alleleFreq'][1][1], subPop[1]['alleleFreq'][2][1],"
               "subPop[1]['alleleFreq'][3][1], subPop[1]['alleleFreq'][4][1], subPop[1]['alleleFreq'][5][1],"
               "subPop[1]['alleleFreq'][6][1], subPop[1]['alleleFreq'][7][1], subPop[1]['alleleFreq'][8][1],"
               "subPop[1]['alleleFreq'][9][1], subPop[1]['alleleFreq'][100][1])", step=1, begin = 73),
       
            #sim.PyOperator(func=lethalEvent)
        ],
        finalOps= sim.Stat(alleleFreq=0, vars=['alleleFreq_sp']),
        gen = 83
    )
    print("HERE")
    pop.vars()

    # print out population size and allele frequency
    for idx, name in enumerate(pop.subPopNames()):
        print('%s (%d): %.4f' % (name, pop.subPopSize(name), 
            pop.dvars(idx).alleleFreq[0][0]))
    
    

    viewVars(pop.vars(), gui=False)
    sim.dump(pop)


    ped = sim.Pedigree(pop);
    print("This is the pedigree stuff")


    # Now sample the individuals
    #sample = drawRandomSample(pop, sizes=[sample_count]*sub_population_count)
    sample = drawRandomSample(pop, sizes=sample_count)

    # Print out the allele frequency data
    sim.stat(sample, alleleFreq=sim.ALL_AVAIL)
    frequencies = sample.dvars().alleleFreq;
    with open(get_filename('freq.txt'), 'w') as freqfile:
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
        
    # =============================================================================
    # MK: Print out mtDNA frequency data in a way that is comparable to the original
    # file.        
    # =============================================================================
    # Print out the mtDNA frequency data
    samplemtDNA_count = 1992 #to grab all the possible final haplotypes; this is limited
    #by te highest number of individuals in a subpopulation, which is 1,992 in this case
    #for subpopulation 1. Once the number of inidividuals is corrected based on the bottleneck,
    #this value here has to be changed.
    samplemtDNA = drawRandomSample(pop, sizes = [samplemtDNA_count]*sub_population_count)
    sim.stat(pop, alleleFreq=sim.ALL_AVAIL)

    last_locus = nb_loci  # Index of the last locus

    with open(get_filename('freq_mtDNA.txt'), 'w') as freqfile:
        index = 0
        for locus in frequencies:
            if locus != last_locus:
                continue

            if len(frequencies[locus]) < 2:
                continue

            print(index, end=' ', file=freqfile)
            index += 1

            for allele in frequencies[locus]:
                print(frequencies[locus][allele], end=' ', file=freqfile)

            print(file=freqfile)

    #The output has values that add up to 1, which is what is expected for a 
    #mtDNA haplotype frequency data.

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
    with open(get_filename('mixfile.txt'), 'w') as mixfile:
        with open(get_filename('haploiso.txt'), 'w') as haplofile:
            print(sub_population_count, nb_loci - nb_ignored_loci, 2, 1, file=mixfile)
            print("sex, haplotype, carbon, nitrogen, native_ground", file=haplofile);
            for i in range(0, nb_loci - nb_ignored_loci):
                print('Loc', i+1, sep='_', file=mixfile);
            for individual in sample.individuals():
                genotype = individual.genotype();
                print(1 if individual.sex() == 1 else 0,
                      genotype[nb_loci],
                      individual.info('carbon'),
                      individual.info('nitrogen'),
                          int(individual.info('native_breeding_ground')),
                      file=haplofile, sep=' ')
                print(int(individual.info('native_breeding_ground')+1), end=' ', file=mixfile)
                for i in range(0, nb_loci):
                    if i not in monoallelic_loci:
                        print(genotype[i]+1, genotype[i+nb_loci+1]+1, ' ', end='', sep='', file=mixfile)
                print(file=mixfile);


if __name__ == '__main__':
    runSimulation(
        [pop_staticPhaseEnd],
        minMatingAge,
        maxMatingAge,
        gen,
        mitochondrial_file = "mtdna_1.single.txt",
        snp_file = 'snp_1.single.txt'
    )