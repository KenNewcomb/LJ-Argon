## System imports
import sys
import math

## Local imports
from classes.simulation import Simulation

# Instanstiate a new Simulation Object, with 864 particles.
sim = Simulation(864)

# Give each atom a random position in 3D-space
sim.randomizePositions()
	
# Apply random velocities to particles.
sim.applyBoltzmannDist()

# Allow them to move according to some time step

# Integrate, find new positions and momenta

# Repeat until finished.

# Analyze results (create velocity autocorrelation function and from this, compute material properties)
