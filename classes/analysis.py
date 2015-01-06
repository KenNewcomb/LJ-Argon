"""analysis.py: A class that provides analysis of the output data."""

class Analysis:
    
    sigma = 3.4e-10 # sigma in Lennard-Jones Potential, meters
    dr = sigma/100 # (1/100)*sigma
    V = (10.229*sigma)**3 # Volume of box
    numAtoms = 864 # Number of atoms
    
    atoms = [] # List of atoms
    nr = [] # n(r), number of particles at radius r
    
    def __init__(self, atoms):
        self.atoms = atoms
    
    def pairDistributionFunction(self):
        """Generates a pair-distribution function on given data"""
        for step in range(0, 5*sigma, dr):
            pass