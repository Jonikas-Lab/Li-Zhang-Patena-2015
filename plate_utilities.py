#! /usr/bin/env python

"""
______
 -- Weronika Patena, 2013
USAGE: _____
"""

# standard library
from __future__ import division
import sys
import unittest
import string
# other packages
# my modules

PLATE_PREFIXES_REMOVE = ['plate']
PLATE_SIZE = 384
ROW_SIZE = 24

def mutant_number_from_plate_well(plate, well):
    """ Get sequential mutant number from plate+well: first plate is 1-384, second is 385-768, etc. """
    # get the plate number - we don't know if it's a string or an int, strip prefixes if string, etc
    if hasattr(plate, 'startswith'):
        for prefix in PLATE_PREFIXES_REMOVE:
            if plate.lower().startswith(prefix):
                plate = plate[len(prefix):]
    plate_number = int(plate)
    # get the well number - adding 1 because the original version is 0-based
    well_row = string.uppercase.index(well[0]) + 1
    well_column = int(well[1:])
    well_number = (well_row-1)*ROW_SIZE + well_column
    # get mutant ID
    return (plate_number-1)*PLATE_SIZE + well_number
    

# TODO should probably merge this with robotic_plate_transfer.Plate_type - the unit-tests are very similar, etc...

class Testing(unittest.TestCase):
    """ Runs unit-tests for this module. """

    def test__mutant_number_from_plate_well(self):
        assert mutant_number_from_plate_well(1, 'A1') == 1
        assert mutant_number_from_plate_well(1, 'A1') == 1
        assert mutant_number_from_plate_well('plate1', 'A1') == 1
        assert mutant_number_from_plate_well('plate1', 'A2') == 2
        assert mutant_number_from_plate_well('plate1', 'A24') == 24
        assert mutant_number_from_plate_well('plate1', 'B1') == 25
        assert mutant_number_from_plate_well('plate1', 'P24') == 384
        assert mutant_number_from_plate_well('plate2', 'A1') == 385
        assert mutant_number_from_plate_well('plate2', 'A2') == 386
        assert mutant_number_from_plate_well('plate3', 'A1') == 769


if __name__=='__main__':
    """ If module is run directly, run tests. """
    print "This is a module for import by other programs - it doesn't do anything on its own.  Running tests..."
    unittest.main()
