"""simulate.py: The main driver for the simulation."""

## System imports
import sys
import math
import os

## Local imports
from classes.simulation import Simulation
from classes.analysis import Analysis
from classes.filewriter import fileWriter
from classes.atom import Atom

# Remove the output file if it exists
try:
    os.remove("output.csv")
except OSError:
    pass

# Number of time steps to run
nSteps = 200

# Instantiate a new simulation object.
sim = Simulation()
analysis = Analysis(sim.getAtoms())
fw = fileWriter()

def printStepInformation(temp):
    print("Current System Temperature: " + str(temp))
    print("-----------------COMPLETED STEP " + str(step+1) + " --------------------")

# Run the simulation for n steps
for step in range(0, nSteps):
    sim.runSimulation(step)
    temp = sim.getTemperature(step)
    printStepInformation(temp)
    
    analysis.updateAtoms(sim.getAtoms())
    analysis.velocityAutocorrelation(step)
analysis.pairDistributionFunction()