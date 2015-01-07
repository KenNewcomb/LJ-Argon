"""simulate.py: The main driver for the simulation."""

## System imports
import sys
import math
import os

## Local imports
from classes.simulation import Simulation
from classes.analysis import Analysis
from classes.filewriter import fileWriter

# Remove the output file if it exists
try:
    os.remove("output.csv")
except OSError:
    pass

# Number of time steps to run
nSteps = 100

# Instantiate a new simulation object.
sim = Simulation()
analysis = Analysis()
fw = fileWriter()

# Run the simulation for n steps
for step in range(0, nSteps):
    sim.runSimulation(step)
    temp = sim.getTemperature
    fw.addTemp(temp)
    analysis.velocityAutocorrelation(temp)
    
    
