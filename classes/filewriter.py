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

    def writeToFile(self, filename, data):
        """Writes data at each step to a file"""
        with open(filename, "a") as output:
            for point in data:
                output.write("%s\n" % point)
