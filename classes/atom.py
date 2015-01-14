"""atom.py: An object that represents a single "atom" to be studied."""

class Atom:
    
    def __init__(self, atomNum):
        # Atom number
        self.num = "Ar"
        
        # Spatial (cartesian) coordinates
        self.x = 0
        self.y = 0
        self.z = 0
        
        # Velocity components
        self.vx = 0
        self.vy = 0
        self.vz = 0
        
        # Force components
        self.fx = 0
        self.fy = 0
        self.fz = 0
