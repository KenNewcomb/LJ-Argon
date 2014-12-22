## System imports
import sys
import math

## Local imports
from classes.simulation import Simulation

### INITIALIZATION ###

# Instanstiate a new Simulation Object, with 864 particles.
sim = Simulation(864)

# Give each atom a random position in 3D-space
sim.randomizePositions()
	
# Apply random velocities to particles.
sim.applyBoltzmannDist()

### MAIN PROGRAM LOOP ###

# Calculate forces on all atoms
sim.updateForces()

# Move the atoms through a time step.
sim.timeStep()

# Analyze results (create velocity autocorrelation function and from this, compute material properties)
