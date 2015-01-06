"""simulation.py: a class that represents a simulation, storing all atomic information and simulation parameters"""

## System imports
import random
import math
import time

## Local imports
from atom import Atom

class Simulation:

    # Simulation parameters/constants 
    kb = 1.380e-23 # Boltzmann (J/K)
    Nav = 6.022e23 # Molecules/mol
    m = (39.95/Nav)*(10**-3) # mass of a single atom , Kg
    e = kb*120 # depth of potential well, J
    sigma = 3.4e-10 # sigma in Lennard-Jones Potential, meters
    rcut = 2.25*sigma # Cutoff radius, meters
    rcutsq = rcut**2 # Square of the cutoff radius.
    T = 90 # Temperature, K
    numAtoms = 864 # Number of atoms to simulate
    lbox = 10.229*sigma # length of the box. (meters)
    dt = 1e-14 # Time step, seconds
    nSteps = 500 # Number of time steps
    realTemp = 0 # System temperature

    atoms = []

    def __init__(self):
        """Creates a simulation with numAtoms"""
        for i in range(0,self.numAtoms):
            self.atoms.append(Atom())

    def assignPositions(self):
        """Places each atom in arbitrary positions in the box."""
        n = 10 # Number of atoms in a direction
        particle = 0 # Particles placed so far
        
        for x in range(0, n):
            for y in range(0, n):
                for z in range(0, n):
                    if (particle < self.numAtoms):
                        self.atoms[particle].x = (x + 0.5) * self.sigma
                        self.atoms[particle].y = (y + 0.5) * self.sigma                 
                        self.atoms[particle].z = (z + 0.5) * self.sigma
                    particle += 1
                    z += 1
                y += 1
            x += 1

    def applyBoltzmannDist(self):
        """Applies Boltzmann distribution to atomic velocities"""
        normDist = []
        scaling_factor = math.sqrt(self.kb*self.T/self.m)

        # Establish Normal Distribution
        for i in range(0, 3*self.numAtoms):
            normDist.append(random.gauss(0,1))
        
        # Apply scaling factor to distribution
        for number in range(0, len(normDist)-1):
            normDist[number] = normDist[number]*scaling_factor
            
        # Distribute velocities
        rand_index = 0
        for atom in range(0, len(self.atoms)-1):
            self.atoms[atom].vx = normDist[rand_index]
            self.atoms[atom].vy = normDist[rand_index+1]
            self.atoms[atom].vz = normDist[rand_index+2]
            rand_index += 3

        for atom in range(0, len(self.atoms)-1):
            self.atoms[atom].xprev = self.atoms[atom].x - self.atoms[atom].vx*self.dt
            self.atoms[atom].yprev = self.atoms[atom].y - self.atoms[atom].vy*self.dt
            self.atoms[atom].zprev = self.atoms[atom].z - self.atoms[atom].vz*self.dt

    def mainLoop(self):
        for step in range(0, self.nSteps):
            self.updateForces()
            self.timeStep()
            self.getTemperature()
            self.resetForces()
            self.scaleTemperature()
            self.writeToFile()
            print("-----------------COMPLETED STEP " + str(step+1) + " --------------------")

    def ljforce(self, atom1, atom2):
        """Calculates the force between two atoms using LJ 12-6 potential"""
        # Calculate distance between two atoms
        dx = self.atoms[atom1].x - self.atoms[atom2].x
        dy = self.atoms[atom1].y - self.atoms[atom2].y
        dz = self.atoms[atom1].z - self.atoms[atom2].z
        
        # Minimum Image Convention
        dx -= self.lbox*round(dx/self.lbox)
        dy -= self.lbox*round(dy/self.lbox)
        dz -= self.lbox*round(dz/self.lbox)
        
        r2 = dx*dx + dy*dy + dz*dz

        if r2 < self.rcutsq:
            fr2 = (self.sigma**2)/r2
            fr6 = fr2**3
            force = 48*self.e*fr6*(fr6 - 0.5)/r2
            
            
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
        for atom in range(0, len(self.atoms)-1):
            # Calculate new positions

            newX = self.atoms[atom].x + self.atoms[atom].vx*self.dt + (self.dt**2)*(self.atoms[atom].fx/self.m)
            newY = self.atoms[atom].y + self.atoms[atom].vy*self.dt + (self.dt**2)*(self.atoms[atom].fy/self.m)
            newZ = self.atoms[atom].z + self.atoms[atom].vz*self.dt + (self.dt**2)*(self.atoms[atom].fz/self.m)

            # Update current velocities
            self.atoms[atom].vx = (newX - self.atoms[atom].x)/(self.dt)
            self.atoms[atom].vy = (newY - self.atoms[atom].y)/(self.dt)
            self.atoms[atom].vz = (newZ - self.atoms[atom].z)/(self.dt)
            
            # Update previous positions
            self.atoms[atom].xprev = self.atoms[atom].x
            self.atoms[atom].yprev = self.atoms[atom].y
            self.atoms[atom].zprev = self.atoms[atom].z
            
            # Update current positions (applying PBC)
            if newX < 0:
                self.atoms[atom].x = newX + self.lbox
                self.atoms[atom].xprev += self.lbox
            elif newX > self.lbox:
                self.atoms[atom].x = newX - self.lbox
                self.atoms[atom].xprev -= self.lbox
            else:
                self.atoms[atom].x = newX
            
            if newY < 0:
                self.atoms[atom].y = newY + self.lbox
                self.atoms[atom].yprev += self.lbox
            elif newY > self.lbox:
                self.atoms[atom].y = newY - self.lbox
                self.atoms[atom].yprev -= self.lbox
            else:
                self.atoms[atom].y = newY
                
            if newZ < 0:
                self.atoms[atom].z = newZ + self.lbox
                self.atoms[atom].zprev += self.lbox
            elif newZ > self.lbox:
                self.atoms[atom].z = newZ - self.lbox
                self.atoms[atom].zprev -= self.lbox 
            else:
                self.atoms[atom].z = newZ

    def resetForces(self):
        """Sets all forces to zero"""
        for atom in range(0, len(self.atoms) -1):
            self.atoms[atom].fx = 0
            self.atoms[atom].fy = 0
            self.atoms[atom].fz = 0
            
    def getTemperature(self):
        """Calculates the current system temperature"""
        sumv2 = 0
        for atom in self.atoms:
            sumv2 += atom.vx**2 + atom.vy**2 + atom.vz**2
        self.realTemp = (self.m/(3*self.numAtoms*self.kb))*sumv2
        print("TEMP: " + str(self.realTemp))
        
    def scaleTemperature(self):
        """Scales the temperature according to desired temperature"""
        if self.realTemp > 100.0 or self.realTemp < 80.0:
            print("Rescaling temperatures...")
            for atom in range(0, len(self.atoms)-1):
                self.atoms[atom].vx *= math.sqrt(self.T/self.realTemp)
                self.atoms[atom].vy *= math.sqrt(self.T/self.realTemp)
                self.atoms[atom].vz *= math.sqrt(self.T/self.realTemp)
                
    def writeToFile(self):
        with open("output.csv", "a") as myfile:
            myfile.write(str(self.realTemp))
            myfile.write("\n")