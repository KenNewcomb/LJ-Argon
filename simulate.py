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

# Number of time steps to run
nSteps = 250

### INITIALIZATION ###

# Instantiate a new simulation object.
sim = Simulation()
# Instantiate an analysis object to analyze the current simulation.
analysis = Analysis(sim.getAtoms())
# Instatiate an object that writes output to the disk.
fw = fileWriter()

### MAIN LOOP ###

# Run the simulation for n steps
for step in range(0, nSteps):
    # Run the simulation for a single step
    sim.runSimulation(step)
    # Add the current state to the analysis object
    analysis.updateAtoms(sim.getAtoms())
    # Calculate the velocity autocorrelation function for this step.
    analysis.velocityAutocorrelation(step)
    
### DATA OUTPUT ###
fw.writeData("temp.csv", sim.getTemperature())
fw.writeData("rdf.csv", analysis.pairDistributionFunction())
fw.writeData("vac.csv", analysis.getVAC())
fw.writeXYZ(sim.atomlist)
