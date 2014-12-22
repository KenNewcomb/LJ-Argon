"""simulation.py: a class that represents a simulation, storing all atomic information and simulation parameters"""

## System imports
import random
import math

## Local imports
from classes import atom

class Simulation:

	# Constants involved in a simulation
	kb = 1.380e-23

	# Simulation/atomic parameters
	m = (39.95*1.6747e-24) # mass of a single atom , Kg
	e = kb*120 # depth of potential well, J
	d = 3.4 # sigma in Lennard-Jones Potential ,angstroms
	rcut = 2.25*d # Cutoff radius , angstroms
	T = 120 # Temperature, K

	# List of atoms
	numAtoms = 0
	atoms = []

	def __init__(self, numAtoms):
		"""Create a simulation with numAtoms"""
		self.numAtoms = numAtoms
		for i in range(0,numAtoms):
			self.atoms.append(atom)
	
	def randomizePositions(self):
		"""Give each atom random positions within the box"""
		for atom in self.atoms:
			atom.x = random.random()*34.78
			atom.y = random.random()*34.78
			atom.z = random.random()*34.78
	
	def applyBoltzmannDist(self):
		"""Apply Boltzmann distribution to atomic velocities"""
		normDist = []
		scaling_factor = math.sqrt(self.kb*self.T/self.m)

		# Establish Normal Distribution
		for i in range(0,3*self.numAtoms):
			normDist.append(random.gauss(0,1))
	
		# Apply scaling factor to distribution
		for number in normDist:
			number = number*scaling_factor
	
		# Distribute velocities	
		rand_index = 0
		for atom in self.atoms:
			atom.vx = normDist[rand_index]
			atom.vy = normDist[rand_index]
			atom.vz = normDist[rand_index]
			rand_index = rand_index + 3
   
       def updateForces(self)
           for atom in self.atoms:
               pass
