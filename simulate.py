import sys
from modules import atomicProperties as ap
from classes import atom
import random

## Basic program structure (to be implemented):

# Create 864 atom objects, store in list.
atoms = []

for i in range(1,864):
	atoms.append(atom)

# Give them random positions and momenta (applying Boltzmann distribution @ T)

for atom in atoms:
	atom.x = random.random()*34.78
	atom.y = random.random()*34.78
	atom.z = random.random()*34.78

# Allow them to move according to some time step


# Integrate, find new positions and momenta

# Repeat until finished.

# Analyze results (create velocity autocorrelation function and from this, compute material properties)
