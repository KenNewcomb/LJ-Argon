"""simulation.py: a class that represents a simulation, storing all atomic information and simulation parameters"""

## System imports
import random
import math

## Local imports
from atom import Atom

class Simulation:

    # Simulation parameters/constants 
    kb = 1.380e-23 # Boltzmann
    m = (39.95*1.6747e-24) # mass of a single atom , Kg
    e = kb*120 # depth of potential well, J
    d = 3.4 # sigma in Lennard-Jones Potential, angstroms
    rcut = 2.25*d # Cutoff radius, angstroms
    rcutsq = rcut**2 # Square of the cutoff radius.
    T = 120 # Temperature, K
    numAtoms = 864 # Number of atoms to simulate
    lbox = 10.229*d # length of the box.
    dt = 10e-14 # Time step, seconds
    
    atoms = []

    def __init__(self):
        """Creates a simulation with numAtoms"""
        for i in range(0,self.numAtoms):
            self.atoms.append(Atom())

    def randomizePositions(self):
        """Gives each atom random positions within the box"""      
        for atom in self.atoms:
            atom.x = random.random()*self.lbox
            atom.y = random.random()*self.lbox
            atom.z = random.random()*self.lbox
        
    def applyBoltzmannDist(self):
        """Applies Boltzmann distribution to atomic velocities"""
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
            atom.vy = normDist[rand_index+1]
            atom.vz = normDist[rand_index+2]
            rand_index = rand_index + 3
            
    def ljforce(self, atom1, atom2):
        """Calculates the force between two atoms using LJ 12-6 potential"""
        # Calculate distance between two atoms
        dx = self.atoms[atom1].x - self.atoms[atom2].x
        dy = self.atoms[atom1].y - self.atoms[atom2].y
        dz = self.atoms[atom1].z - self.atoms[atom2].z
        r2 = dx*dx + dy*dy + dz*dz
        
        # Minimum image convention 
        dx = dx - self.lbox*round(dx/self.lbox)
        dy = dy - self.lbox*round(dy/self.lbox)
        dz = dz - self.lbox*round(dz/self.lbox)
        
        # If atoms are within cutoff radius
        if r2 < self.rcutsq:
            recipr = 1/r2
            recipr6 = r2**3
            force = 48*self.e*recipr*recipr6*(recipr-0.5)
            forcex = force*dx
            forcey = force*dy
            forcez = force*dz
            
            self.atoms[atom1].fx += forcex
            self.atoms[atom2].fx -= forcex
            self.atoms[atom1].fy += forcey
            self.atoms[atom2].fy -= forcey
            self.atoms[atom1].fz += forcez
            self.atoms[atom2].fz -= forcez
        
    def updateForces(self):
        """Calculates the net potential on each atom, applying a cutoff radius"""
        for atom1 in range(0, len(self.atoms)-1):
            for atom2 in range(atom1+1, len(self.atoms)):
                self.ljforce(atom1, atom2)

    def timeStep(self):
        """Moves the system through a given time step, according to the energies"""
        for atom in self.atoms:
            # Calculate new positions
            newX = 2*atom.x - atom.xprev + self.dt**2*atom.fx
            newY = 2*atom.y - atom.yprev + self.dt**2*atom.fy
            newZ = 2*atom.z - atom.zprev + self.dt**2*atom.fz

            # Update current velocities
            atom.vx = (newX - atom.x)/(2*self.dt)
            atom.vy = (newY - atom.y)/(2*self.dt)
            atom.vz = (newZ - atom.z)/(2*self.dt)
            
            # Update previous positions
            atom.xprev = atom.x
            atom.yprev = atom.y
            atom.zprev = atom.z
            
            # Update current positions (applying PBC)
            if newX < 0:
                atom.x = newX + self.lbox
            elif newX > self.lbox:
                atom.x = newX - self.lbox
            else:
                atom.x = newX
            
            if newY < 0:
                atom.y = newY + self.lbox
            elif newY > self.lbox:
                atom.y = newY - self.lbox
            else:
                atom.y = newY
                
            if atom.z < 0:
                atom.z = newZ + self.lbox
            elif newZ > self.lbox:
                atom.z = newZ - self.lbox
            else:
                atom.z = newZ
                
        print("Atom 1 x position: " + str(self.atoms[0].x))
        print("Atom 1 x force: " + str(self.atoms[0].fx))