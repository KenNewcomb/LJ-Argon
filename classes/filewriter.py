# filewriter.py: A module writes the output of a simluation to a file.

import os

class fileWriter:
    
    def __init__(self):
        #Remove the output file if it exists
        try:
            os.remove("output.csv")
        except OSError:
            pass
    
    def addTemp(self, temp):
        """Adds a temperature to the current list of temperatures"""
        self.temperatures.append(temp)

    def writeTemperatures(self, temperatures):
        with open("output.csv", "a") as output:
            for temperature in temperatures:
                output.write("%s\n" % temperature)

