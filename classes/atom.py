# An object that represents a single "atom" to be studied.
# Each atom has Cartesian coordinates x, y, and z, as well as their respective
# components of velocity (vx, vy, and vz).

class Atom:
    
    def __init__(self):
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