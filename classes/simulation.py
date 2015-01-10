"""simulation.py: a class that represents a simulation, storing all atomic information and simulation parameters"""

## System imports
import random
import math
import copy

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
    temp = 90 # Temperature, K
    numAtoms = 864 # Number of atoms to simulate
    lbox = 10.229*sigma # length of the box. (meters)
    dt = 1e-14 # Time step, seconds
    nSteps = 10 # Number of time steps
    currentTemp = 0 # System temperature
    epot = 0 # Potential energy

    atoms = []
    temperatures = []
    
    def __init__(self):
        """Creates a simulation with numAtoms"""
        for i in range(0,self.numAtoms):
            self.atoms.append(Atom())
        self.assignPositions()
        self.applyBoltzmannDist()
        self.correctMomenta()
        
    def assignPositions(self):
        """Places each atom in arbitrary positions in the box."""
        n = int(math.ceil(self.numAtoms**(1.0/3.0))) # Number of atoms in a direction
        particle = 0 # Particles placed so far
        
        for x in range(0, n-1):
            for y in range(0, n):
                for z in range(0, n):
                    if (particle < self.numAtoms):
                        self.atoms[particle].x = x * self.sigma*1.1
                        self.atoms[particle].y = y * self.sigma             
                        self.atoms[particle].z = z * self.sigma
                        particle += 1

    def applyBoltzmannDist(self):
        """Applies Boltzmann distribution to atomic velocities"""
        normDist = []
        scaling_factor = math.sqrt(self.kb*self.temp/self.m)

        # Establish Normal Distribution
        for i in range(0, 3*self.numAtoms):
            normDist.append(random.gauss(0,1))
      
        # Apply scaling factor to distribution
        for number in range(0, 3*self.numAtoms):
            normDist[number] = normDist[number]*scaling_factor
            
        # Distribute velocities
        for atom in range(0, self.numAtoms):
            self.atoms[atom].vx = normDist[atom*3]
            self.atoms[atom].vy = normDist[atom*3+1]
            self.atoms[atom].vz = normDist[atom*3+2]

    def correctMomenta(self):
        sumvx = 0
        sumvy = 0
        sumvz = 0
        
        for atom in range(0, self.numAtoms):
            sumvx += self.atoms[atom].vx
            sumvy += self.atoms[atom].vy
            sumvz += self.atoms[atom].vz
        
        for atom in range(0, self.numAtoms):
            self.atoms[atom].vx -= sumvx/self.numAtoms
            self.atoms[atom].vy -= sumvy/self.numAtoms
            self.atoms[atom].vz -= sumvz/self.numAtoms
            
        
    def runSimulation(self, step):
        self.updateForces()
        self.verletIntegration()
        self.updateTemperature()
        self.resetForces()
        print("Current System Temperature: " + str(self.currentTemp))
        print("-----------------COMPLETED STEP " + str(step+1) + " --------------------")
        # After 100 steps, scale the temperature by a factor of (Tdesired/T(t))^1/2
        if step > 100:
            self.scaleTemperature()
        
    def updateForces(self):
        """Calculates the net potential on each atom, applying a cutoff radius"""
        for atom1 in range(0, self.numAtoms-1):
            for atom2 in range(atom1+1, self.numAtoms):
                self.calculateForce(atom1, atom2)
                
        # Multiply by constants        
        for atom in range(0, self.numAtoms):
            self.atoms[atom].fx *= 48*self.e
            self.atoms[atom].fy *= 48*self.e
            self.atoms[atom].fz *= 48*self.e
            
    def calculateForce(self, atom1, atom2):
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
            force = fr6*(fr6 - 0.5)/r2
            
            forcex = force*dx
            forcey = force*dy
            forcez = force*dz
            
            self.atoms[atom1].fx += forcex
            self.atoms[atom2].fx -= forcex
            self.atoms[atom1].fy += forcey
            self.atoms[atom2].fy -= forcey
            self.atoms[atom1].fz += forcez
            self.atoms[atom2].fz -= forcez
            
    def verletIntegration(self):
        """Moves the system through a given time step, according to the energies"""
        for atom in range(0, self.numAtoms):
            
            # Update velocities
            self.atoms[atom].vx += (self.atoms[atom].fx/self.m)*self.dt
            self.atoms[atom].vy += (self.atoms[atom].fy/self.m)*self.dt
            self.atoms[atom].vz += (self.atoms[atom].fz/self.m)*self.dt

            # Update positions
            newX = self.atoms[atom].x + self.atoms[atom].vx*self.dt
            newY = self.atoms[atom].y + self.atoms[atom].vy*self.dt
            newZ = self.atoms[atom].z + self.atoms[atom].vz*self.dt

            # Update current positions (applying PBC)
            if newX < 0:
                self.atoms[atom].x = newX + self.lbox
            elif newX > self.lbox:
                self.atoms[atom].x = newX - self.lbox
            else:
                self.atoms[atom].x = newX
            
            if newY < 0:
                self.atoms[atom].y = newY + self.lbox
            elif newY > self.lbox:
                self.atoms[atom].y = newY - self.lbox
            else:
                self.atoms[atom].y = newY
                
            if newZ < 0:
                self.atoms[atom].z = newZ + self.lbox
            elif newZ > self.lbox:
                self.atoms[atom].z = newZ - self.lbox
            else:
                self.atoms[atom].z = newZ

    def resetForces(self):
        """Sets all forces to zero"""
        for atom in range(0, self.numAtoms):
            self.atoms[atom].fx = 0
            self.atoms[atom].fy = 0
            self.atoms[atom].fz = 0
            
    def updateTemperature(self):
        """Calculates the current system temperature"""
        sumv2 = 0
        for atom in self.atoms:
            sumv2 += atom.vx**2 + atom.vy**2 + atom.vz**2
        self.currentTemp = (self.m/(3*self.numAtoms*self.kb))*sumv2
        self.temperatures.append(self.currentTemp)
    
    def getAtoms(self):
        return copy.deepcopy(self.atoms)
    
    def getTemperature(self):
        return copy.deepcopy(self.temperatures)
        
    def scaleTemperature(self):
        """Scales the temperature according to desired temperature"""
        if self.currentTemp > 100.0 or self.currentTemp < 80.0:
            print("Rescaling velocities...")
            for atom in range(0, self.numAtoms):
                self.atoms[atom].vx *= math.sqrt(self.temp/self.currentTemp)
                self.atoms[atom].vy *= math.sqrt(self.temp/self.currentTemp)
                self.atoms[atom].vz *= math.sqrt(self.temp/self.currentTemp)
