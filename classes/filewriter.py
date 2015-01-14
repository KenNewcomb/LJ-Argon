# filewriter.py: A module writes the output of a simluation to a file.

import os

class fileWriter:
    
    def __init__(self):
        #Remove the output file if it exists
        try:
            os.remove("temperatures.csv")
        except OSError:
            pass
        
        try:
            os.remove("rdf.csv")
        except OSError:
            pass

        try:
            os.remove("vac.csv")
        except OSError:
            pass
        
        try:
            os.remove("ljargon.xyz")
        except OSError:
            pass

    def writeData(self, filename, data):
        """Writes data at each step to a file"""
        with open(filename, "a") as output:
            for point in data:
                output.write("%s\n" % point)

    def writeXYZ(self, xyzs):
        """Writes positional data to a .xyz file"""
        with open("ljargon.xyz", "a") as output:
            output.write("864\n") # Number of atoms
            output.write("Ar\n") # Name of Molecule
            for frame in xyzs:
                for atom in frame:
                    output.write("%s %s %s %s\n" % (atom.num, atom.x, atom.y, atom.z))
