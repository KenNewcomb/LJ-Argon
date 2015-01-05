## System imports
import sys
import math

## Local imports
from classes.simulation import Simulation

### INITIALIZATION ###

# Instanstiate a new simulation object.
sim = Simulation()

# Initially place each atom on a simple cubic lattice
sim.assignPositions()
	
# Apply random velocities to particles.
sim.applyBoltzmannDist()

### MAIN PROGRAM LOOP ###
sim.mainLoop()