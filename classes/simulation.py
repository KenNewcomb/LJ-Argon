"""simulation.py: a class that represents a simulation, storing all atomic information and simulation parameters"""

## System imports
import random
import math

## Local imports
from atom import Atom

class Simulation:

    # Constants involved in a simulation
    kb = 1.380e-23 # Boltzmann

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
        """Creates a simulation with numAtoms"""
        self.numAtoms = numAtoms
        for i in range(0,numAtoms):
            self.atoms.append(Atom())

    def randomizePositions(self):
        """Gives each atom random positions within the box"""      
        for i in range(0, self.numAtoms):
            self.atoms[i].x = random.random()*34.78
            self.atoms[i].y = random.random()*34.78
            self.atoms[i].z = random.random()*34.78
        
    def applyBoltzmannDist(self):
        """Applies Boltzmann distribution to atomic velocities"""
        normDist = []
        scaling_factor = math.sqrt(self.kb*self.T/self.m)
        print("scaling factor " + str(scaling_factor))

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
            atom.vy = normDist[rand_index+1]
            atom.vz = normDist[rand_index+2]
            rand_index = rand_index + 3
            
        print(self.atoms[0].vx)
        print(self.atoms[1].vx)
            
    def ljforce(self, atom1, atom2):
        """Calculates the force between two atoms using LJ 12-6 potential"""
        dx = atom1.x - atom2.x
        dy = atom1.y - atom2.y
        dz = atom1.z - atom2.z
        r2 = dx*dx + dy*dy + dz*dz
        #print("r2: " + str(r2))
        
    def updateForces(self):
        """Calculates the net potential on each atom, applying a cutoff radius"""
        for atom1 in range(0, len(self.atoms)-1):
            for atom2 in range(atom1+1, len(self.atoms)):
                self.ljforce(self.atoms[atom1], self.atoms[atom2])