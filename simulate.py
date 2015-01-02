## System imports
import sys
import math

## Local imports
from classes.simulation import Simulation

### INITIALIZATION ###

# Instanstiate a new simulation object.
sim = Simulation()

# Give each atom a random position in 3D-space
sim.randomizePositions()
	
# Apply random velocities to particles.
sim.applyBoltzmannDist()

### MAIN PROGRAM LOOP ###

for steps in range(0,100):
    # Calculate forces on all atoms
    sim.updateForces()
    # Move the atoms through a time step.
    sim.timeStep()

# Analyze results (create velocity autocorrelation function and from this, compute material properties)