"""analysis.py: A class that provides analysis of the output data."""

import math

class Analysis:
    
    sigma = 3.4e-10 # sigma in Lennard-Jones Potential, meters
    dr = sigma/10 # (1/100)*sigma
    V = (10.229*sigma)**3 # Volume of box
    numAtoms = 864 # Number of atoms
    
    originalAtoms = [] # List of atoms
    currentAtoms = []
    nr = [] # n(r), number of particles at radius r
    velacfinit = 0 # Velocity autocorrlation function at time = 0
    velacf = 0 # velocity autocorrelation function at a time step
    lbox = 10.229*sigma # Length of box
    
    velaclist = [] # The velocity autocorrelation function for a simulation
    
    def __init__(self, atoms):
        """Initalizize the analysis with the atoms in their initial state"""
        self.originalAtoms = atoms
        
    def updateAtoms(self, atoms):
        """Update the current state of the atoms"""
        self.currentAtoms = atoms
        
    def pairDistributionFunction(self):
        """Generates a pair-distribution function on given data"""
        atom_counts = []
        atom_count = 0
        cur_r = 0
        
        for step in range(0, 50):
            for atom in range(1, self.numAtoms):
                dx = self.currentAtoms[0].x - self.currentAtoms[atom].x
                dy = self.currentAtoms[0].y - self.currentAtoms[atom].y
                dz = self.currentAtoms[0].z - self.currentAtoms[atom].z

                # Minimum Image Convention
                dx -= self.lbox*round(dx/self.lbox)
                dy -= self.lbox*round(dy/self.lbox)
                dz -= self.lbox*round(dz/self.lbox)
                
                r2 = dx*dx + dy*dy + dz*dz
                r = math.sqrt(r2)
                
                # If the atom is within the range
                if (r < (cur_r + self.dr)) and (r > cur_r):
                    atom_count += 1

            # Increment by dr    
            cur_r += self.dr
            
            # Add the counts for that radius
            atom_counts.append((self.V*atom_count/self.numAtoms)/(4*math.pi*r2*self.dr))
            
            # Reset atom count for next radius
            atom_count = 0
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
        