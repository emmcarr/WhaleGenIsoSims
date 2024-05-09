#!/usr/bin/env python3
# coding=utf-8
"""check simulations"""

import csv
import math
from pathlib import Path

#age,ind_id,father_id,mother_id,nitrogen,carbon,feeding_ground,native_breeding_ground,migrate_to,sex,aff
#33.0,18835.0,0.0,0.0,6.063689437535627,14.072635396536256,0.0,0.0,0.0,M,U

# set up carbon and nitrogen values and bounds as 3*sd
mean_C = [16.7, 20.5, 26.0] 
variance_C = [3.24, 2.89, 3.0]

mean_N = [5.5, 8.7, 5.0]
variance_N = [0.25, 0.49, 0.4]

CARBON, NITROGEN = {}, {}
for i in (0, 1, 2):
    CARBON[i] = (
        mean_C[i] - 3*math.sqrt(variance_C[i]),  # lower bound
        mean_C[i] + 3*math.sqrt(variance_C[i])   # upper bound
    )
    NITROGEN[i] = (
        mean_N[i] - 3*math.sqrt(variance_N[i]),  # lower bound
        mean_N[i] + 3*math.sqrt(variance_N[i])   # upper bound
    )



def read(filename):
    # n.b. the simulator spits out all haplo data is _1 and _2, so this will
    # lose all but the (last) column due to the name clashes. These will be the SNP haplos
    # so that should be ok.
    with open(filename, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for i in reader:
            # do some type casting and things here to make life easier
            # handle ints
            for key in ('ind_id', 'father_id', 'mother_id'):
                i[key] = None if i[key] == '0.0' else int(i[key].replace(".0", ""))
            
            for key in ('feeding_ground', 'native_breeding_ground', 'migrate_to'):
                i[key] = int(i[key].replace(".0", ""))
            
            # handle floats
            for key in ('nitrogen', 'carbon', 'age'):
                i[key] = float(i[key])
            yield i



def check_sex(i):
    """Check that sex is M or F"""
    if i['sex'] not in ('M', 'F'):
        print("ERROR: Unusual sex obtained for individual %d: %s" % (i['ind_id'], i['sex']))


def check_feeding_ground(i):
    """Check feeding_ground"""
    if i['feeding_ground'] not in (0, 1, 2):
        print("ERROR: Unusual feeding_ground obtained for individual %d: %s" % (i['ind_id'], i['feeding_ground']))


def check_native_breeding_ground(i):
    """Check native_breeding_ground"""
    if i['native_breeding_ground'] != 0:
        print("ERROR: Unusual native_breeding_ground obtained for individual %d: %s" % (i['ind_id'], i['native_breeding_ground']))


def check_nitrogen(i):
    """
    Check that nitrogen is not too wild (within 4 s.d. of mean)

    mean_N = [5.5, 8.7, 5.0]
    variance_N = [0.25, 0.49, 0.4]
    """
    expected_range = NITROGEN[i['feeding_ground']]
    if not expected_range[0] < i['nitrogen'] < expected_range[1]:
        print("ERROR: Nitrogen outside bounds for individual %d: %0.4f, ground=%d" % (i['ind_id'], i['nitrogen'], i['feeding_ground']))


def check_carbon(i):
    """
    Check that carbon is not too wild (within 4 s.d. of mean)
    
    mean_C = [16.7, 20.5, 26.0] 
    variance_C = [3.24, 2.89, 3.0]
    """
    expected_range = CARBON[i['feeding_ground']]
    if not expected_range[0] < i['carbon'] < expected_range[1]:
        print("ERROR: Carbon outside bounds for individual %d: %0.4f, ground=%d" % (i['ind_id'], i['carbon'], i['feeding_ground']))


def check_parental_haplos(i, population):
    parental1 = [population.get(i['mother_id'], {}).get('_1'), population.get(i['father_id'], {}).get('_1')]
    parental2 = [population.get(i['mother_id'], {}).get('_2'), population.get(i['father_id'], {}).get('_2')]
    
    parental1 = [h for h in parental1 if h]
    parental2 = [h for h in parental2 if h]

    if parental1 and i['_1'] not in parental1:
        print("ERROR: haplogroup mismatch %d: %r != %r" % (i['ind_id'], i['_1'], parental1))
    
    if parental2 and i['_2'] not in parental2:
        print("ERROR: haplogroup mismatch %d: %r != %r" % (i['ind_id'], i['_2'], parental2))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Does something.')
    parser.add_argument("filename", help='filename', type=Path)
    args = parser.parse_args()
    
    population = {indiv['ind_id']: indiv for indiv in read(args.filename)}

    for indiv_id, indiv in population.items():
        check_sex(indiv)
        check_feeding_ground(indiv)
        check_native_breeding_ground(indiv)
        check_nitrogen(indiv)
        check_carbon(indiv)
        check_parental_haplos(indiv, population)