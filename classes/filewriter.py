# filewriter.py: A module writes the output of a simluation to a file.

class fileWriter:

    temperatures = []
    
    def addTemp(self, temp):
        """Adds a temperature to the current list of temperatures"""
        self.temperatures.append(temp)

    def writeToFile(self):
        with open("output.csv", "a") as output:
            for temperature in temperatures:
                output.write("%s\n" % temperature)

