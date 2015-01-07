import unittest
from classes.simulation import Simulation

class Test_Simulation(unittest.TestCase):
    
    def test_generate_particles(self):
        """Tests for the creation of numAtoms in a simulation"""
        sim = Simulation()
        self.assertEqual(sim.numAtoms, 864)
        

if __name__ == '__main__':
    unittest.main()