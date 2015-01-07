"""simulate.py: The main driver for the simulation."""

## System imports
import sys
import math
import os

## Local imports
from classes.simulation import Simulation
from classes.analysis import Analysis

### INITIALIZATION ###
try:
    os.remove("output.csv")
except OSError:
    pass
# Instantiate a new simulation object.
sim = Simulation()

# Run the simulation.
sim.runSimulation()
