#!/usr/bin/env python3
# coding=utf-8
"""check simulations"""

import csv
from pathlib import Path

from check_simulation import read

def pp(gen, i):
    #print(i)
    return "%2d\t%-6d\t%1s\t%2d\t%6d\t%6d\t%0.3f\t%0.3f\t%d\t%d\t%d\t%s\t%s" % (
        gen, i['ind_id'], i['sex'], i['age'], 
        # set mother and father ids to -1 if the individual is from generation 0
        i['mother_id'] if i['mother_id'] else -1,
        i['father_id'] if i['father_id'] else -1,
        i['carbon'], i['nitrogen'],
        i['feeding_ground'], i['native_breeding_ground'], i['migrate_to'],
        i['_1'], i['_2']
    ) 

def tostr(i):
    return "<%-6d: %s %dy (C=%0.3f, N=%0.3f, FG=%s, BG=%s, M=%s, H=%s,%s)>" % (
        i['ind_id'], i['sex'], i['age'], i['carbon'], i['nitrogen'],
        i['feeding_ground'], i['native_breeding_ground'],
        i['migrate_to'],
        i['_1'], i['_2']
    )
    
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Does something.')
    parser.add_argument("simulation_id", help='simulationid')
    parser.add_argument("individual", help='individual', type=int)
    args = parser.parse_args()
    
    known_children = set()
    mother, father = None, None
    print("Gen\tID\tSex\tAge\tMother\tFather\tC\tN\tFG\tBG\tM\tH1\tH2")
    for i in range(1, 83):
        filename = Path("results") / ("%s_dump_gen_%d.csv" % (args.simulation_id, i))
        population = {indiv['ind_id']: indiv for indiv in read(filename)}
        if args.individual in population:
            
            # find parents if we haven't done so yet
            if not mother:
                mother = population.get(population[args.individual]['mother_id'])
                if mother:
                    print('<- MOTHER', tostr(mother))

            if not father:
                father = population.get(population[args.individual]['father_id'])
                if father:
                    print('<- FATHER', tostr(father))
            
            
            print(pp(i, population[args.individual]))
            # find offspring
            parentkey = 'father_id' if population[args.individual]['sex'] == 'M' else 'mother_id'
            offspring = [
                o for o in population.values() if o[parentkey] == args.individual
                if o['ind_id'] not in known_children
            ]
            for o in offspring:
                print('-> OFFSPRING', tostr(o))
                known_children.add(o['ind_id'])
            
        else:
            print("%2d" % i)
    

