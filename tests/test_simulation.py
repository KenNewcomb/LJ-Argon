import unittest
from classes.simulation import Simulation

class Test_Simulation(unittest.TestCase):
    
    sim = Simulation()
    
    def test_generate_particles(self):
        """Tests for the creation of numAtoms in a simulation"""
        self.assertEqual(len(self.sim.atoms), 864)
        
    def test_number_of_particles_placed_on_lattice(self):
        """Test to ensure no two atoms are on the same lattice positions"""
        self.sim.assignPositions()
        coordlist = []
        
        for atom in range(0, self.sim.numAtoms):
            coord = (self.sim.atoms[atom].x, self.sim.atoms[atom].y, self.sim.atoms[atom].z)
            coordlist.append(coord)
            
        self.assertEqual(len(coordlist), len(set(coordlist)))
if __name__ == '__main__':
    unittest.main()