## System imports
import sys
import random
import math

## Local imports
from modules import atomicProperties as ap
from classes import atom

# Create 864 atom objects, store in list.
num_particles = 864
atoms = []
for i in range(0,864):
	atoms.append(atom)

# Give each atom a random position in 3D-space
for atom in atoms:
	atom.x = random.random()*34.78
	atom.y = random.random()*34.78
	atom.z = random.random()*34.78
	
# Generate 3N (N=864) random numbers from normal distribution.
norm_dist = []
for i in range(0,3*num_particles):
	norm_dist.append(random.gauss(0,1))

# Apply scaling factor (sqrt(kT/m))
scaling_factor = math.sqrt(ap.kb*ap.T/ap.m)
for number in norm_dist:
	number = number*scaling_factor

# Apply random velocities to particles.
rand_index = 0
for atom in atoms:
	atom.vx = norm_dist[rand_index]
	atom.vy = norm_dist[rand_index+1]
	atom.vz = norm_dist[rand_index+2]
	rand_index = rand_index+3

# Allow them to move according to some time step


# Integrate, find new positions and momenta

# Repeat until finished.

# Analyze results (create velocity autocorrelation function and from this, compute material properties)
