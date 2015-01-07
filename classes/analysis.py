"""analysis.py: A class that provides analysis of the output data."""

class Analysis:
    
    sigma = 3.4e-10 # sigma in Lennard-Jones Potential, meters
    dr = sigma/100 # (1/100)*sigma
    V = (10.229*sigma)**3 # Volume of box
    numAtoms = 864 # Number of atoms
    
    originalAtoms = [] # List of atoms
    currentAtoms = []
    nr = [] # n(r), number of particles at radius r
    velacfinit = 0
    velacf = 0 # velocity autocorrelation function
    
    def __init__(self, atoms):
        """Initalizize the analysis with the atoms in their initial state"""
        self.originalAtoms = atoms
        
    def updateAtoms(self, atoms):
        """Update the current state of the atoms"""
        self.currentAtoms = atoms
        
    def pairDistributionFunction(self):
        """Generates a pair-distribution function on given data"""
        atom_counts = []
        for step in range(0, 5*sigma, dr):
            for atom in range(0, self.numAtoms):
                dx = self.atoms[atom1].x - self.atoms[atom2].x
                dy = self.atoms[atom1].y - self.atoms[atom2].y
                dz = self.atoms[atom1].z - self.atoms[atom2].z
        
                # Minimum Image Convention
                dx -= self.lbox*round(dx/self.lbox)
                dy -= self.lbox*round(dy/self.lbox)
                dz -= self.lbox*round(dz/self.lbox)
        
                r2 = dx*dx + dy*dy + dz*dz
    
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
        else:   
            for atom in range(0, self.numAtoms):
                vx += self.originalAtoms[atom].vx * self.currentAtoms[atom].vx
                vy += self.originalAtoms[atom].vy * self.currentAtoms[atom].vy
                vz += self.originalAtoms[atom].vz * self.currentAtoms[atom].vz
            self.velacf += vx + vy + vz
            self.velacf /= self.numAtoms*self.velacfinit
            print(self.velacf)
            
        self.velacf = 0
        