# -*- coding: utf-8 -*-
"""
Created on Fri Nov 24 13:28:45 2023

@author: marin
"""

# =============================================================================
# Down here: the code merged with forwardTrajectory.py that worked
# =============================================================================

pop.evolve(
    initOps = [sim.InitSex(), sim.IdTagger(), sim.PyOperator(func=init_native_breeding_grounds)] +
                   [sim.InitGenotype(subPops = sub_population_names[i], freq=haplotype_frequencies[i], loci=[nb_loci]) for i in range(0, sub_population_count)] + 
                   #this means that the mtDNA is the very last locus
                   [sim.InitGenotype(subPops = sub_population_names[i], freq=[snp[n][i], 1-snp[n][i]], loci=[n]) for i in range(0, sub_population_count) for n in range(0, nb_loci-1)],
    # increase age by 1
    preOps = [sim.InfoExec('age += 1')],
    matingScheme=sim.ControlledRandomMating(
        ops=[sim.Recombinator(rates=0.01)],
        loci=5, alleles=1, freqFunc=traj.func()),
    postOps = [
        sim.Stat(structure=range(1), subPops=sub_population_names, suffix='_AB', step=10),
        sim.PyEval(r"'Fst=%.3f \n' % (F_st_AB)", step=10)
    ],
    gen = 80
)



traj = simulateForwardTrajectory(N=[2000, 4000], fitness=None,
    beginGen=0, endGen=100, beginFreq=[0.2, 0.3],
    endFreq=[[0.1, 0.11], [0.2, 0.21]])
# 
#traj.plot('log/forwardTrajectory.png', set_ylim_top=0.5,
#    plot_c_sp=['r', 'b'], set_title_label='Simulated Trajectory (forward-time)')
pop = sim.Population(size=[2000, 4000], loci=10, infoFields='fitness')
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitGenotype(freq=[0.8, 0.2], subPops=0),
        sim.InitGenotype(freq=[0.7, 0.3], subPops=1),
        sim.PyOutput('Sp0: loc2\tloc5\tSp1: loc2\tloc5\n'),
    ],
    matingScheme=sim.ControlledRandomMating(
        ops=[sim.Recombinator(rates=0.01)],
        loci=5, alleles=1, freqFunc=traj.func()),
    postOps=[
        sim.Stat(alleleFreq=[2, 5], vars=['alleleFreq_sp'], step=20),
        sim.PyEval(r"'%.2f\t%.2f\t%.2f\t%.2f\n' % (subPop[0]['alleleFreq'][2][1],"
            "subPop[0]['alleleFreq'][5][1], subPop[1]['alleleFreq'][2][1],"
            "subPop[1]['alleleFreq'][5][1])", step=20),
        
    ],
    gen = 101
)


# =============================================================================
# Code that replaced randomMating for ControlledrandomMating, now obsolete
# =============================================================================
pop.evolve(
    initOps = [sim.InitSex(), sim.IdTagger(), sim.PyOperator(func=init_native_breeding_grounds)] +
                   [sim.InitGenotype(subPops = sub_population_names[i], 
                                     freq=haplotype_frequencies[i], 
                                     loci=[nb_loci]) for i in range(0, sub_population_count)] +
                   [sim.InitGenotype(subPops = sub_population_names[i], 
                                     freq=[snp[n][i], 1-snp[n][i]], 
                                     loci=[n]) for i in range(0, sub_population_count) for n in range(0, nb_loci-1)],
    # increase age by 1
    
    preOps = [sim.InfoExec('age += 1')],
    #matingScheme = sim.HeteroMating([
    matingScheme = sim.ControlledRandomMating(
        #MK: removing options until it stops breaking:
            
        # age <= maxAge, copy to the next generation (weight=-1)
        # subPops is a list of tuples that will participate in mating. The tuple is a pair (subPopulation, virtualSubPopulation)
        # First, we propagate (clone) all individuals in all subpopulations (and all VSPs except the ones who are now in the VSP of deceased individuals) to the next generation
        # sim.CloneMating(ops=[sim.CloneGenoTransmitter(chroms=[0,1])],
          #                  subPops=[(sub_population, 6) for sub_population in range(0, sub_population_count)],
           #                 weight=-1),
        
        # Then we simulate random mating only in VSP 1 (ie reproductively mature individuals) within subpopulation (breeding/winter grounds)
       # sim.RandomMating(ops=[sim.MitochondrialGenoTransmitter(),
        #                          sim.MendelianGenoTransmitter(),
         #                         sim.IdTagger(),
          #                        sim.InheritTagger(mode=sim.MATERNAL, infoFields=['feeding_ground']),
           #                       sim.InheritTagger(mode=sim.MATERNAL, infoFields=['native_breeding_ground']),
            #                      sim.PedigreeTagger()],
             #                subPops=[(sub_population, 7) for sub_population in range(0, sub_population_count)],
              #               weight=1)],
            #subPopSize=configure_new_population_size,
               #                               ),

                              subPopSize=configure_new_population_size,
                              freqFunc=traj.func()),
    
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
        sim.PyEval(r"'Fst=%.3f \n' % (F_st_AB)", step=10)
    ],
    gen = 80
)

sim.dump(pop)

#sim.dump(pop, width=3, loci=[], subPops=[(sim.ALL_AVAIL, sim.ALL_AVAIL)], max=1000, structure=False);
#return



# =============================================================================
# Does it work with the heteromating, trajectory, bottleneck, a.k.a., 
# merging them all together?
# =============================================================================
# Emma's initial code, once again. It works!
pop.evolve(
    initOps = [sim.InitSex(), sim.IdTagger(), sim.PyOperator(func=init_native_breeding_grounds)] +
                   [sim.InitGenotype(subPops = sub_population_names[i], freq=haplotype_frequencies[i], loci=[nb_loci]) for i in range(0, sub_population_count)] +
                   [sim.InitGenotype(subPops = sub_population_names[i], freq=[snp[n][i], 1-snp[n][i]], loci=[n]) for i in range(0, sub_population_count) for n in range(0, nb_loci-1)],
    # increase age by 1
    preOps = [sim.InfoExec('age += 1')],
    matingScheme = sim.HeteroMating([
        # age <= maxAge, copy to the next generation (weight=-1)
        # subPops is a list of tuples that will participate in mating. The tuple is a pair (subPopulation, virtualSubPopulation)
        # First, we propagate (clone) all individuals in all subpopulations (and all VSPs except the ones who are now in the VSP of deceased individuals) to the next generation
        sim.CloneMating(ops=[sim.CloneGenoTransmitter(chroms=[0,1])],
                            subPops=[(sub_population, 6) for sub_population in range(0, sub_population_count)], ##6 because it's dividing the pop into the VSPs based on age and removing
                            #the dead ones
                            weight=-1), 
        # Then we simulate random mating only in VSP 1 (ie reproductively mature individuals) within subpopulation (breeding/winter grounds)
        sim.RandomMating(subPopSize=model,
            ops=[sim.MitochondrialGenoTransmitter(),
                                  sim.MendelianGenoTransmitter(),
                                  sim.IdTagger(),
                                  sim.InheritTagger(mode=sim.MATERNAL, infoFields=['feeding_ground']),
                                  sim.InheritTagger(mode=sim.MATERNAL, infoFields=['native_breeding_ground']),
                                  sim.PedigreeTagger()],
                             subPops=[(sub_population, 7) for sub_population in range(0, sub_population_count)],
                             weight=1),
        sim.ControlledRandomMating(freqFunc=traj.func(),
                                   weight=1)] 
        ),
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
        sim.PyEval(r"'Fst=%.3f \n' % (F_st_AB)", step=10), #Print Fst every 10 steps
        sim.Stat(alleleFreq=[2, 5], vars=['alleleFreq_sp'], step=10), #added this now, to
        #calculate the allele frequencies in selected loci
        sim.PyEval(r"'%.2f\t%.2f\t%.2f\t%.2f\n' % (subPop[0]['alleleFreq'][2][1],"
            "subPop[0]['alleleFreq'][5][1], subPop[1]['alleleFreq'][2][1],"
            "subPop[1]['alleleFreq'][5][1])", step=10) #added this now, to
            #print out the allele frequencies in selected loci
       
    ],
    finalOps= sim.Stat(alleleFreq=0, vars=['alleleFreq_sp']),
    gen = model4.num_gens
    #80
)


# =============================================================================
# Printing only mitochondrial frequencies:
# =============================================================================

pop.evolve(
    initOps = [sim.InitSex(), sim.IdTagger(), sim.PyOperator(func=init_native_breeding_grounds)] +
                   [sim.InitGenotype(subPops = sub_population_names[i], freq=haplotype_frequencies[i], loci=[nb_loci]) for i in range(0, sub_population_count)] +
                   [sim.InitGenotype(subPops = sub_population_names[i], freq=[snp[n][i], 1-snp[n][i]], loci=[n]) for i in range(0, sub_population_count) for n in range(0, nb_loci-1)],
    # increase age by 1
    preOps = [sim.InfoExec('age += 1')],
    matingScheme = sim.HeteroMating([
        # age <= maxAge, copy to the next generation (weight=-1)
        # subPops is a list of tuples that will participate in mating. The tuple is a pair (subPopulation, virtualSubPopulation)
        # First, we propagate (clone) all individuals in all subpopulations (and all VSPs except the ones who are now in the VSP of deceased individuals) to the next generation
        sim.CloneMating(ops=[sim.CloneGenoTransmitter(chroms=[0,1])],
                            subPops=[(sub_population, 6) for sub_population in range(0, sub_population_count)],
                            weight=-1), #if weights are negative, they are multiplied to their parental subpopulation;
            #EX: parental pop: (500, 1000), weight -2, next generation: (1000, 2000). For weight -1, it keeps the number of individuals
            #from the parental generation.
            #ALSO: if there is a mix of negative and positive weights, the negative will be processed first.
        # Then we simulate random mating only in VSP 1 (ie reproductively mature individuals) within subpopulation (breeding/winter grounds)
        sim.RandomMating(subPopSize=model4,
            ops=[sim.MitochondrialGenoTransmitter(),
                                  sim.MendelianGenoTransmitter(),
                                  sim.IdTagger(),
                                  sim.InheritTagger(mode=sim.MATERNAL, infoFields=['feeding_ground']),
                                  sim.InheritTagger(mode=sim.MATERNAL, infoFields=['native_breeding_ground']),
                                  sim.PedigreeTagger()],
                             subPops=[(sub_population, 7) for sub_population in range(0, sub_population_count)],
                             weight=1),
        sim.ControlledRandomMating(freqFunc=traj.func(),
                                   weight=1)] #we decided to keep the same weight as the
        #mitochondrial transmitter.
        ),
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
        sim.Stat(numOfMales=True, 
                 haploFreq=[], 
                 begin = 70, step = 75, end = 85),     
        sim.PyEval(r"'Fst=%.3f \n' % (F_st_AB)", step=10), #Print Fst every 10 steps
        sim.Stat(alleleFreq=[1, 2, 3, 4, 5, 6, 7, 8, 9, 100], vars=['alleleFreq_sp'], step=10), #added this now, to
        #calculate the allele frequencies in selected loci
       sim.PyEval(r"'%.2f\t%.2f\n' % ("
           "subPop[0]['alleleFreq'][100][1]," 
           "subPop[1]['alleleFreq'][100][1])", step=10),
      # sim.PyEval(r"'%.2f' % haploFreq[0]", step=10)

        #to print out the allele frequencies in selected loci
       
    ],
    finalOps= sim.Stat(alleleFreq=0, vars=['alleleFreq_sp']),
    gen = model4.num_gens
    )