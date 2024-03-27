# -*- coding: utf-8 -*-
"""
Created on Fri Nov 24 13:28:01 2023

@author: marin
"""

# =============================================================================
# Model 1: a simple (and ugly) model
# =============================================================================
model1 = demo.MultiStageModel([
    demo.InstantChangeModel(T=108, #number of generations
        # start with an ancestral population of size 30000
        N0=(30000, 'Ancestral'),
        # change population size at 105 and 107
        G=[105, 107], 
        # change to population size 40 and back to 1000
        NG=[(40, 'Bottleneck'), (1000, 'Post-Bottleneck')]),
    demo.ExponentialGrowthModel(
        T=3, #since they started recovering from the bottleneck 
        # split the population into two subpopulations
        N0=[(500, 'A'), (500, 'B')],
        # expand to size 1000 and 1000
        NT=[1000, 1000])]
    )

model1.init_size #returns the initial population size
model1.info_fields
model1.num_gens #sum of T=108 and T=3

# get a visual presentation of the demographic model
model1.plot('log/simplemodel.png',
    title='Model 1: A simple (ugly) model')

# =============================================================================
# Working on a more realistic model
# =============================================================================
#Useful link: https://simupop.sourceforge.net/manual_release/build/userGuide_ch7_sec3.html?highlight=multistagemodel

model2 = demo.MultiStageModel([
    #initialize the population with 60 individuals so I can use the haplotype 
    #frequencies file to initialize their haplotypes
    demo.ExponentialGrowthModel(
        T=10, #10 initial generations ("burn in" similar to Emma's)
        N0=[(60, 'A'), (60, 'B')], #splits into two subpopulations, one keeps 100% of the individuals
        #and the second one has 60
        NT=[15000, 15000] #changes the population to 15000 each in 100 generations
        ),
    demo.AdmixtureModel(model = ('CGF', #continuous gene flow
                            0, 1, 0.9), T=100), #mixes 10% of individuals from 
    #population 1 to population 2 for 10 generations
    demo.InstantChangeModel(T=108, #number of generations
        # start with an ancestral population of size 30000
        #N0=(30000, 'Ancestral'), #do not need to specify because it follows the previous stage
        # change population size at 105 and 107
        G=[105, 107], 
        # change to population size 40 and back to 1000
        NG=[(40, 'Bottleneck'), (1000, 'Post-Bottleneck')]),
    demo.ExponentialGrowthModel(
        T=3, #since they started recovering from the bottleneck 
        # split the population into two subpopulations
        N0=[(500, 'C'), (500, 'D')],
        # expand to size 1000 and 1000
        NT=[1000, 1000])]
    )

model2.init_size #returns the initial population size
model2.info_fields
model2.num_gens #sum of T=108 and T=3

# get a visual presentation of the demographic model
model2.plot('log/averagemodel.png',
    title='Model 2: An average model')

# =============================================================================
# Same as before, but with more layers
# =============================================================================
# To allow trajectory to be calculated, it has to be a number between 80 and 1000,
# from my tries. Here, I will start at 70:
    # Gen 70: 1830-1855 maximum amount of whales
    # Gen 75: 1930-1955 lowest amount (bottleneck)
    # Gen 78: 1980-2005 recovering 
    # Gen 85: 2010-2030 current(ish) generation

model3 = demo.MultiStageModel([
        demo.ExponentialGrowthModel(
        T=70, #1 generation
        N0=(20000, 'Ancestral'), 
        NT=[30000] 
        ),
    #AdmixtureModel(model = ('CGF', #continuous gene flow
     #                       0, 1, 0.9), T=100), #mixes 10% of individuals from 
    #population 1 to population 2 for 10 generations
    demo.InstantChangeModel(T=5, # generation on year 1855
                       G=[3, 4], #generarion on year 1855, year 1905
    NG=[(3000, "Start of bottleneck"), (200)]),        
    demo.InstantChangeModel(T=8, # generation on year 1930
        # start with an ancestral population of size 30000
        #N0=(30000, 'Ancestral'), #do not need to specify because it follows the previous stage
        # change population size at 105 and 107
        G=[6, 7], #generarion on year 1955, year 1980
        # change to population size 40 and back to 1000
        NG=[(150), (1000, 'Post-Bottleneck')]),
    demo.ExponentialGrowthModel(
        T=2, #since they started recovering from the bottleneck 
        # split the population into two subpopulations
        N0=[(500), (500)],
        # expand to size 1000 and 1000
        NT=[1000, 1000])]
    )

# get a visual presentation of the demographic model
model3.plot('log/demoModel.png',
    title='Model 3: A simple model with more datapoints')

model3.init_size #returns the initial population size
model3.info_fields
model3.num_gens 

# =============================================================================
# Now that it all works, is there a way to create a more linear whaling curve?
# =============================================================================

model4 = demo.MultiStageModel([
        demo.ExponentialGrowthModel(
        T=75, #75 generations
        N0=(25000, 'Ancestral'), 
        NT=[30000] 
        ),
        demo.ExponentialGrowthModel(
        T=5, #5 generations
        N0=(30000, 'Start of bottleneck'), 
        NT=[3000] 
        ),
    #AdmixtureModel(model = ('CGF', #continuous gene flow
     #                       0, 1, 0.9), T=100), #mixes 10% of individuals from 
    #population 1 to population 2 for 10 generations
   #Changing the instant change to a exponential change:
    demo.ExponentialGrowthModel(
    T=1, #1 generation
    N0=(3000, 'Bottleneck'), 
    NT=[200] 
    ), 
   demo.ExponentialGrowthModel(
   T=1, #1 generation
   N0=(200, 'Bottleneck'), 
   NT=[150] 
   ),
  demo.ExponentialGrowthModel(
  T=1, #1 generation
  N0=(150, 'Bottleneck'), 
  NT=[40] 
  ),
   demo.ExponentialGrowthModel(
   T=1, #1 generation
   N0=(40, 'Bottleneck'), 
   NT=[150] 
   ),
   demo.ExponentialGrowthModel(
   T=1, #1 generation
   N0=(150, 'Recovery'), 
   NT=[1000] 
   ),
   demo.ExponentialGrowthModel(
   T=2, #2 generations
   N0=(1000, 'Recovery'), 
   NT=[2500] 
   )])

    # get a visual presentation of the demographic model
model4.plot('log/demoModel.png',
    title='Model 4: A simple model with more exponential changes')

model4.init_size #returns the initial population size
model4.info_fields
model4.num_gens 



# =============================================================================
# Model 5: model 4, but with two subpopulations
# =============================================================================
model5 = demo.MultiStageModel([
        demo.ExponentialGrowthModel(
        T=75, #65 generations
        N0=[12500, 12500], 
        NT=[15000, 15000] 
        ),
        demo.ExponentialGrowthModel(
        T=5, #5 generations
        N0=[(1500, 'A'), (1500, 'B')] , 
        NT=[(1500, 'A'), (1500, 'B')] 
        ),
    demo.AdmixtureModel(model = ('CGF', #continuous gene flow
                           0, 1, 0.9), T=10), #mixes 10% of individuals from 
    #population 1 to population 2 for 10 generations
   #Changing the instant change to a exponential change:
    demo.ExponentialGrowthModel(
    T=1, #1 generation
    N0=[(1500, 'A'), (1500, 'B')] , 
    NT=[(100, 'A'), (100, 'B')]  
    ),
   demo.ExponentialGrowthModel(
   T=1, #1 generation
   N0=[(100, 'A'), (100, 'B')],
   NT=[(75, 'A'), (75, 'B')]  
   ),
  demo.ExponentialGrowthModel(
  T=1, #1 generation
  N0=[(75, 'A'), (75, 'B')], 
  NT=[(20, 'A'), (20, 'B')] 
  ),
   demo.ExponentialGrowthModel(
   T=1, #1 generation
   N0=[(20, 'A'), (20, 'B')], 
   NT=[(75, 'A'), (75, 'B')]
   ),
   demo.ExponentialGrowthModel(
   T=1, #1 generation
   N0=[(75, 'A'), (75, 'B')], 
   NT=[(500, 'A'), (500, 'B')]
   ),
   demo.ExponentialGrowthModel(
   T=2, #2 generations
   N0=[(500, 'A'), (500, 'B')], 
   NT=[(1250, 'A'), (1250, 'B')]
   )])
    
    # get a visual presentation of the demographic model
model5.plot('log/demoModel.png',
    title='Model 5: A simple model with more exponential changes and subpopulations')

model5.init_size #returns the initial population size
model5.info_fields
model5.num_gens 


# =============================================================================
# Model 6: model based on events
# =============================================================================


#Using bottleneck.py,
#written by Diane Bailleul, with the help of Bo Peng and Solenn Stoeckel.
from simuPOP.demography import *



model6 = EventBasedModel(
    N0=([15000]*2),
    events=[ResizeEvent(at=76, sizes=[1500, 1500]),
            ResizeEvent(at=77, sizes=[100, 100]),
            ResizeEvent(at=78, sizes=[75, 75]),
            ResizeEvent(at=79, sizes=[20, 20]),
            ResizeEvent(at=80, sizes=[75, 75]),
            ResizeEvent(at=81, sizes=[500, 500]),
            ResizeEvent(at=82, sizes=[1069, 1069])
        ]
)

#This event is not able to be visualized on a plot and it also gives a num_gens of 
# -1, for some reason. I checked the original code (by Diane Bailleul, mentioned 
#above) and theirs also returns a -1 generation, so I think it should be ok.

#model6.plot('log/demoModel.png',
 #           title = 'Model 6: A model with events')

model6.init_size #returns the initial population size
model6.info_fields
model6.num_gens 


# =============================================================================
# Model 7: a bottleneck based on simuBottlenecks_v2.py, kindly shared by
#  Juhana Kammonen
# =============================================================================



#Gen 4000: 1830-1855: ~ 30,000 individuals
#Gen 4001: 1855-1880: ~ 3,000
#Gen 4002: 1880-1905: ~ 200
#Gen 4003: 1905-1930: ~ 150
#Gen 4004: 1930-1955: ~ 40
#Gen 4005: 1955-1980: ~ 150
#Gen 4006: 1980-2005: ~ 1000 
#Gen 4007: 2005-2030: ~ 2,139 (2009)
#Gen 4008: 2010-2030: ~ 2500 (maybe, at least?)



# GENERATION OPTIONS (1 gen = 25 years)

staticPhaseEnd = 4000 # 9000 BP - begin steady growth at latest
firstDecline	= 4001 # 5750 BP - decline begins at population peak
firstDeclineContinues = 4002 # the decline evens out before bottleneck
firstDeclineContinues2 = 4003 # the decline evens out before bottleneck
Bottleneck	= 4004 # 4100 BP - first bottleneck begins  
firstRecovery	= 4005 # 3800 BP - first bottleneck ends
secondRecovery	= 4006 # population growth passes a thousand
highestsinceBottleneck = 4007 # self-explanatory
current = 4008 # 2010 - 2030 estimates

firstCheck      = 4000 # population peak
secondCheck     = 4001 # when bottleneck starts happening
thirdCheck      = 4004 # lowest number - bottleneck
fourthCheck    = 4006 # population growth passes a thousand
fifthCheck     = 4007 # population just before the current one
sixthCheck     = 4008 # current generation until 2030s

# =============================================================================
# MK Model 8
# All the other deleted models are in Removed_Models_Whale.py
# =============================================================================

#Gen 1000: 1830-1855: ~ 30,000 individuals
#Gen 1001: 1855-1880: ~ 3,000
#Gen 1002: 1880-1905: ~ 200
#Gen 1003: 1905-1930: ~ 150
#Gen 1004: 1930-1955: ~ 40
#Gen 1005: 1955-1980: ~ 150
#Gen 1006: 1980-2005: ~ 1000 
#Gen 1007: 2005-2030: ~ 2,139 (2009)
#Gen 1008: 2010-2030: ~ 2500 (maybe, at least?)

#Gen 75: 1830-1855: ~ 30,000 individuals
#Gen 76: 1855-1880: ~ 3,000
#Gen 77: 1880-1905: ~ 200
#Gen 78: 1905-1930: ~ 150
#Gen 79: 1930-1955: ~ 40
#Gen 80: 1955-1980: ~ 150
#Gen 81: 1980-2005: ~ 1000 
#Gen 82: 2005-2030: ~ 2,139 (2009)
#Gen 83: 2010-2030: ~ 2500 (maybe, at least?)


model8 = demo.MultiStageModel([
        demo.LinearGrowthModel(
            #up until staticPhaseEnd
        T=75,
        N0=[10000, 10000], 
        NT=[15000, 15000] 
        ),
       #Changing the instant change to a exponential change:
    demo.ExponentialGrowthModel(
    T=1, #first decline
    NT=[(1500), (1500)]  
    ),
   demo.ExponentialGrowthModel(
   T=1, #1 generation
   NT=[(100), (100)]  
   ),
  demo.ExponentialGrowthModel(
  T=1, #1 generation
   NT=[(75), (75)] 
  ),
   demo.ExponentialGrowthModel(
   T=1, #1 generation
  # N0=[(75), (75)], 
   NT=[(20), (20)]
   ),
   demo.ExponentialGrowthModel(
   T=1, #1 generation
   #N0=[(20), (20)], 
   NT=[(75), (75)]
   ),
   demo.ExponentialGrowthModel(
   T=1, #1 generation
   #N0=[(75), (75)], 
   NT=[(500), (500)]
   ),
   demo.LinearGrowthModel(
   T=2, #2 generations
  # N0=[(500), (500)], 
   NT=[(1250), (1250)]
   )])
    
    # get a visual presentation of the demographic model
model8.plot('log/demoModel.png',
    title='Model 8: A simple model with more exponential changes and subpopulations')

model8.init_size #returns the initial population size
model8.info_fields
model8.num_gens 


# =============================================================================
# Model 9: trying to make sure the popSize changes over pop.evolve
# =============================================================================
#From the vignette:
#An instant population growth model that evolves a population from size N0 to NT 
#for T generations with population size changes at generation G to NT.
#If G is a list, multiple population size changes are allowed. In that case, a list 
#(or a nested list) of population size should be provided to parameter NT. 
#Both N0 and NT supports fixed (an integer), dynamic (keep passed poulation size) 
#and proportional (an float number) population size. Optionally, one or more operators 
#(e.g. a migrator) ops can be applied to population. Required information fields by 
#these operators should be passed to parameter infoFields. If removeEmpty option is set to True, 
#empty subpopulation will be removed. This option can be used to remove subpopulations.


model9 = demo.InstantChangeModel(T=83, #generations
                                 N0=[(15000, 15000)], #initial size
                                 G=[76, 77, 78, 79, 80, 81, 82], #generations where the size changes
                                 NG=[[(1500), (1500)],[(100), (100)], [(75), (75)],
                                     [(20), (20)], [(75), (75)], [(500), (500)],
                                     [(1050), (1050)]],
                                 #NG=[3000, 200, 150, 40, 150, 1000, 2139], 
                                 ops=[], 
                                 infoFields=[], 
                                 removeEmptySubPops=False)

model9.plot('log/demoModel.png',
    title='Model 9: A simple model with instant changes')

model9.init_size #returns the initial population size
model9.info_fields
model9.num_gens 


# =============================================================================
# Model 10: Inspired by bottleneck.py (Diane Bailleul)
# =============================================================================
model10 = EventBasedModel(
    N0=([15000]*2), 
    events=[
        ResizeEvent(at=76, sizes=1500),
        ResizeEvent(at=77, sizes=100),
        ResizeEvent(at=78, sizes=75),
        ResizeEvent(at=79, sizes=20),
        ResizeEvent(at=80, sizes=75),
        ResizeEvent(at=81, sizes=500),
        ResizeEvent(at=82, sizes=1050)
        ]
    )


model10.init_size #returns the initial population size
model10.num_gens 



# =============================================================================
# Editing Juhanna's code
# =============================================================================

firstCheck = staticPhaseEnd
secondCheck = firstDecline
thirdCheck = secondDecline
fourthCheck = thirdDecline
fifthCheck = bottleNeck
sixthCheck = firstRecovery
seventhCheck = secondRecovery
eigthCheck = currentGen
ninethCheck = hopefulFuture



    
    
def removeOverspill(pop):

    gen = pop.dvars().gen
    
    if gen <= staticPhaseEnd:
        gen += 1 # increment gen
        nextSubPopSizes = [eval(argList[9]), eval(argList[10])] # (SAA, MUU) growth until max
        maxAllowed = sum(nextSubPopSizes)
        numRemovables = pop.popSize() - maxAllowed 
        overspillLethalEvent(pop, numRemovables)

    elif gen <= firstBalance:
        gen += 1 # increment gen
        nextSubPopSizes = [eval(argList[11]), eval(argList[12])] # (SAA, MUU) growth until max
        maxAllowed = sum(nextSubPopSizes)
        numRemovables = pop.popSize() - maxAllowed 
        overspillLethalEvent(pop, numRemovables)

    
    elif gen < firstDecline:
        gen += 1
        nextSubPopSizes = [eval(argList[13]), eval(argList[14]), eval(argList[15])] #(SAA, SW, NE)
        maxAllowed = sum(nextSubPopSizes)
        numRemovables = pop.popSize() - maxAllowed
        overspillLethalEvent(pop, numRemovables)
        
    elif gen < firstMinimum:
        gen += 1
        nextSubPopSizes = [eval(argList[16]), eval(argList[17]), eval(argList[18])]
        maxAllowed = sum(nextSubPopSizes)
        numRemovables = pop.popSize() - maxAllowed
        overspillLethalEvent(pop, numRemovables)
        
    elif gen < secondBalance:
        gen += 1
        nextSubPopSizes = [eval(argList[19]), eval(argList[20]), eval(argList[21])]
        maxAllowed = sum(nextSubPopSizes)
        numRemovables = pop.popSize() - maxAllowed
        overspillLethalEvent(pop, numRemovables)
        
    elif gen < firstRecovery:
        gen += 1
        nextSubPopSizes = [eval(argList[22]), eval(argList[23]), eval(argList[24])]
        maxAllowed = sum(nextSubPopSizes)
        numRemovables = pop.popSize() - maxAllowed
        overspillLethalEvent(pop, numRemovables)
        
    elif gen < thirdBalance:
        gen += 1
        nextSubPopSizes = [eval(argList[25]), eval(argList[26]), eval(argList[27])]
        maxAllowed = sum(nextSubPopSizes)
        numRemovables = pop.popSize() - maxAllowed
        overspillLethalEvent(pop, numRemovables)
        
    elif gen < secondDecline:
        gen += 1
        nextSubPopSizes = [eval(argList[28]), eval(argList[29]), eval(argList[30])]
        maxAllowed = sum(nextSubPopSizes)
        numRemovables = pop.popSize() - maxAllowed
        overspillLethalEvent(pop, numRemovables)
        
    else: #final rise - growth from bottleneck until present
        gen += 1
        nextSubPopSizes = [eval(argList[31]), eval(argList[32]), eval(argList[33])]
        maxAllowed = sum(nextSubPopSizes)
        numRemovables = pop.popSize() - maxAllowed
        overspillLethalEvent(pop, numRemovables)
        
    return True



def lethalEvent(pop): # the idea here is to affect a portion of the population with
                      # a some kind of an lethal effect that "kills" individuals by simply removing them.
                      # This is a step closer to the "more realistic" simulation model.

    gen = pop.dvars().gen

    #Could be something like this:
    indices = range(pop.popSize()) # a list of individual indices
    numAffected = int(.15 * pop.popSize()) # this would make the affection ratio 15% of current population

    affected = random.sample(indices, numAffected) # picks numAffected indices

    #finally remove selected individuals
    pop.removeIndividuals(affected)

    return True


pop = Population(size=[250, 250],
                 ploidy=2,
                 loci=[631, 0, 16], # 631 loci in mitochondrial DNA, one for each base - 16 loci in Y-chromosome
                 chromTypes=[CUSTOMIZED, CHROMOSOME_X, CHROMOSOME_Y],
                 infoFields=['age', 'group', 'migrate_to'], ancGen=0



# Evolve the POPULATION (self-evolving populations are a simuPOP 1.0.0+ feature)
# ----------------------
pop.evolve(

    initOps = [PyOperator(func=initPop),
               PyOperator(func=initGroups),
               PyOperator(func=sample),
               PyOperator(func=printVsp)],
    
    preOps = [
	# split into two subpops at generation 401 (7000 BP)
	SplitSubPops(subPops=[1], proportions=[0.333,1-0.333], at=401),
        PyOperator(func=loadBgPops),
        InfoExec('age += 1', begin=2),
        PyOperator(func=cleanUp, begin=2),

        # Migration from neighbouring populations in the 2-subpop-phase:
        PyOperator(insertMigration, param=[0], begin=0, end=400),
        PyOperator(insertMigration, param=[1], begin=0, end=400),

        # Migration from neighbouring populations in the 3-subpop-phase:
        PyOperator(insertMigration, param=[0], begin=401),
        PyOperator(insertMigration, param=[1], begin=401),
        PyOperator(insertMigration, param=[2], begin=401),

        # Internal migration (2-subpop-phase):
        Migrator(rate=[
            [0, malnear], # males from subpop 0 to subpops 0 and 1
            [0, femnear], # females from subpop 0 to subpops 0 and 1
            [malnear, 0], # males from subpop 1 to subpops 1 and 0
            [femnear, 0]  # females from subpop 1 to subpops 1 and 0
            ],
            mode=BY_PROPORTION,
            subPops=[(0,4), (0,5), (1,4), (1,5)],
            end=400
            ),

        # Internal migration (3-subpop-phase):
        Migrator(rate=[
            [0, malnear, malfar], # males from subpop 0 to subpops [0, 1, 2]
            [0, femnear, femfar], # females from subpop 0 to subpops [0, 1, 2]
            [malnear, 0, malnear], # ...
            [femnear, 0, femnear],  # ...
            [malfar, malnear, 0],
            [femfar, femnear, 0], # females from subpop 2 to subpops [0, 1, 2]
            ],
            mode=BY_PROPORTION,
            subPops=[(0,4), (0,5), (1,4), (1,5), (2,4), (2,5)],
            begin=401
            ),

        PyOperator(func=removeOverspill), # randomly reduce population size so that parent generation fits.
        PyOperator(func=printNumSubPops),
        
        # mutation call for mitochondrial DNA sequences. Takes deletions into account too:
        MatrixMutator(rate=[
            [0, u/4, (u/4)*kappa, u/4, 0],
            [u/4, 0, u/4, (u/4)*kappa, 0],
            [(u/4)*kappa, u/4, 0, u/4, 0],
            [u/4, (u/4)*kappa, u/4, 0, 0],
            [0, 0, 0, 0, 1]
            ],
                      loci = range(631)),
        
        # mutation call for Y_chromosomal microsatellites
        StepwiseMutator(rates=float(argList[8]), loci=range(631,647)), # mutation rates equal
        
        # debug printing:
        #PyEval(r"'Generation: %d\n' %(gen) "),
        PyOperator(func=printSize),

        PyOperator(func=sample, at=[firstCheck, secondCheck, thirdCheck, fourthCheck, fifthCheck, sixthCheck, seventhCheck, eighthCheck, ninthCheck]),
        PyOperator(func=dumpSequences, at=ninthCheck),
        ],

    # mating scheme similar to that of Savonian expansion simulations:
    matingScheme = HeteroMating(matingSchemes=[
                ## BO: I see why you have to use setMatingScheme, which is not necessary when you use
                ## the new ALL_AVAIL feature in 1.0.1
                ##
                ## Juhana: setMatingScheme has been removed from script

                 # age <= maxAge, copy to the next generation (weight=-1).
                 # no need to worry about too old individuals, they're removed before mating.
                 CloneMating(subPops=[(ALL_AVAIL, 0), (ALL_AVAIL, 1), (ALL_AVAIL, 2)],
                             ops=[
                             ## BO: because CloneGenoTransmitter does not
                             ## by default handle customized chromosomes, an
                             ## explicit list of chromosomes have to be used.
                             CloneGenoTransmitter(chroms=[0,2])
                             ],
                ## BO: I do not see weight = -1, this is important because other wise the number of cloned
                ## individuals will be determined by relative number of individuals compared to randomMating guys.
                             weight=-1),

                 # random mating for individuals in mating ages (vsp 1 in every subpop)
                 RandomMating(subPops=[(ALL_AVAIL, 1)],
                              ops=[MitochondrialGenoTransmitter(),
                                   MendelianGenoTransmitter(),
                                   #IdTagger(), PedigreeTagger(), ParentsTagger()],
                                   ],
                              numOffspring=(POISSON_DISTRIBUTION, 2))
                 ]
                , subPopSize=demo),

    postOps=[
	TicToc(at=[0,1099], output=">>>serial_running_times.txt"),
        PyOperator(func=lethalEvent),
        Stat(popSize=True, subPops=[(ALL_AVAIL, 0), (ALL_AVAIL, 1), (ALL_AVAIL, 2), (ALL_AVAIL, 3)]),
        #PyEval(r"'Size of age groups: %s\n' % (','.join(['%d' % x for x in subPopSize]))")
        ],
    
    finalOps = [PyOutput('Simulation complete\n'), Dumper(max=10, structure=True)],
    
    gen = numGens
    )
