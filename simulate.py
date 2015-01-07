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

# Initially place each atom on a simple cubic lattice
sim.assignPositions()
	
# Apply random velocities to particles.
sim.applyBoltzmannDist()

### MAIN PROGRAM LOOP ###
sim.mainLoop()
