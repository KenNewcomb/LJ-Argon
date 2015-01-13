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

    def writeTemperatures(self, temperatures):
        """Writes the system temperatures at each step to a file"""
        with open("temperatures.csv", "a") as output:
            for temperature in temperatures:
                output.write("%s\n" % temperature)
    
    def writeRDF(self, rdfs):
        """Writes the radial distribution function to a file"""
        with open("rdf.csv", "a") as output:
            for rdf in rdfs:
                output.write("%s\n" % rdf)
    
    def writeVAC(self, vacs):
        """Writes the velocity autocorrelation function to a file"""
        with open("vac.csv", "a") as output:
            for vac in vacs:
                output.write("%s\n" % vac)
    
    