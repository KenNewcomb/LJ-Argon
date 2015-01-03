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

for steps in range(0,100):
    # Calculate forces on all atoms
    sim.updateForces()
    # Move the atoms through a time step.
    sim.timeStep()
    sim.getTemperature()

# Analyze results (create velocity autocorrelation function and from this, compute material properties)