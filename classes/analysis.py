"""analysis.py: A class that provides analysis of the output data."""

import math
import numpy
import matplotlib.pyplot as plot

class Analysis:
    
    kb = 1.380e-23 # Boltzmann (J/K)
    sigma = 3.4e-10 # sigma in Lennard-Jones Potential, meters
    dr = sigma/10 # (1/100)*sigma
    V = (10.229*sigma)**3 # Volume of box
    numAtoms = 864 # Number of atoms
    dt = 1e-14
    originalAtoms = [] # List of atoms
    currentAtoms = []
    nr = [] # n(r), number of particles at radius r
    velacfinit = 0 # Velocity autocorrlation function at time = 0
    velacf = 0 # velocity autocorrelation function at a time step
    lbox = 10.229*sigma # Length of box
    
    velaclist = [] # The velocity autocorrelation function for a simulation
    radiuslist = []
    timelist = []
    
    def __init__(self, atoms):
        """Initalizize the analysis with the atoms in their initial state"""
        self.originalAtoms = atoms
        
    def updateAtoms(self, atoms):
        """Update the current state of the atoms"""
        self.currentAtoms = atoms
        
    def pairDistributionFunction(self):
        """Generates a pair-distribution function on given data"""
        atom_counts = [0]*50
        cur_r = 0
        
        print("Generating rdf..."),
        
        for atom1 in range(0, self.numAtoms-1):
            for atom2 in range(1, self.numAtoms):
                dx = self.currentAtoms[atom1].x - self.currentAtoms[atom2].x
                dy = self.currentAtoms[atom1].y - self.currentAtoms[atom2].y
                dz = self.currentAtoms[atom1].z - self.currentAtoms[atom2].z

                # Minimum Image Convention
                dx -= self.lbox*round(dx/self.lbox)
                dy -= self.lbox*round(dy/self.lbox)
                dz -= self.lbox*round(dz/self.lbox)
                
                r2 = dx*dx + dy*dy + dz*dz
                r = math.sqrt(r2)
                
                # If the atom is within the 
                for radius in range(0, 50):
                    if (r < ((radius+1)*self.dr)) and (r > radius*self.dr):
                        atom_counts[radius] += 1
            
        # Multiply by constants
        for radius in range(1, 50):
            atom_counts[radius] *= (self.V/self.numAtoms**2)/(4*math.pi*((radius*self.dr)**2)*self.dr)
        print("done.")    
        return(atom_counts)        
                    
    def velocityAutocorrelation(self, step):
        vx = 0
        vy = 0
        vz = 0
        if step == 0:
            for atom in range(0, self.numAtoms):
                vx += self.originalAtoms[atom].vx * self.currentAtoms[atom].vx
                vy += self.originalAtoms[atom].vy * self.currentAtoms[atom].vy
                vz += self.originalAtoms[atom].vz * self.currentAtoms[atom].vz
            self.velacfinit += vx + vy + vz
            self.velacfinit /= self.numAtoms
            self.velaclist.append(self.velacfinit)
        else:   
            for atom in range(0, self.numAtoms):
                vx += self.originalAtoms[atom].vx * self.currentAtoms[atom].vx
                vy += self.originalAtoms[atom].vy * self.currentAtoms[atom].vy
                vz += self.originalAtoms[atom].vz * self.currentAtoms[atom].vz
            self.velacf += vx + vy + vz
            self.velacf /= self.numAtoms*self.velacfinit
            self.velaclist.append(self.velacf)
            self.velacf = 0
    
    def getVAC(self):
        return self.velaclist

    def plotRDF(self):
        """Plots the RDF on the screen"""
        pass
        rdf = numpy.loadtxt("rdf.csv")
        for radius in range(0,50):
            self.radiuslist.append(radius*self.dr)
        plot.figure()
        plot.plot(self.radiuslist,rdf)
        plot.show()
        

    def plotVAC(self, nSteps):
       """Plots the VAC on the screen"""
       pass
       vac = numpy.loadtxt("vac.csv")
       vac[0] = 1
       for time in range(0, nSteps):
           self.timelist.append(float(time) * self.dt)
       plot.figure()
       plot.plot(self.timelist, vac)
       plot.show()

    def plotEnergy(self, temperatures, potentials, nSteps):
        """Plots the kinetic, potential, and total energies of the system"""
        # Temperatures --> Energies
        KE = []
        for temp in temperatures:
            KE.append(3*self.numAtoms*self.kb*temp/2)
        
        # Generate a list of steps
        steplist = []
        for time in range(0, nSteps):
            steplist.append(float(time))
        
        # Generate a total energy function
        etot = []
        for energy in range(0, nSteps):
            etot.append(KE[energy] + potentials[energy])
            
        plot.figure()
        plot.plot(steplist, KE, steplist, potentials, steplist, etot)
        plot.show()