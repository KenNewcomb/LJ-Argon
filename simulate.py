## System imports
import sys
import math
import os

## Local imports
from classes.simulation import Simulation

### INITIALIZATION ###
try:
    os.remove("output.csv")
except OSError:
    pass
# Instanstiate a new simulation object.
sim = Simulation()

# Initially place each atom on a simple cubic lattice
sim.assignPositions()
	
# Apply random velocities to particles.
sim.applyBoltzmannDist()

### MAIN PROGRAM LOOP ###
sim.mainLoop()