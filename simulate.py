"""simulate.py: The main driver for the simulation."""

from classes.simulation import Simulation
from classes.analysis import Analysis
from classes.filewriter import fileWriter
import sys

# Number of time steps to run
if len(sys.argv) != 2:
    print("Usage: python simulate.py nSteps")
    exit(0)
else:
    nSteps = int(sys.argv[1])

### INITIALIZATION ###

# Instantiate a new simulation object.
sim = Simulation()
# Instantiate an analysis object to analyze the current simulation.
analysis = Analysis(sim.getAtoms())
# Instatiate an object that writes output to the disk.
fw = fileWriter()

### MAIN LOOP ###

# Run the simulation for n steps
for step in range(0, nSteps):
    # Run the simulation for a single step
    sim.runSimulation(step)
    # Add the current state to the analysis object
    analysis.updateAtoms(sim.getAtoms())
    # Calculate the velocity autocorrelation function for this step.
    analysis.velocityAutocorrelation(step)
    # Write the positional data to a file
    fw.writeXYZ(sim.getAtoms())
    
### DATA OUTPUT ###
fw.writeData("temp.csv", sim.temperatures)
fw.writeData("rdf.csv", analysis.pairDistributionFunction())
fw.writeData("vac.csv", analysis.getVAC())

analysis.plotRDF()
analysis.plotVAC(nSteps)
analysis.plotEnergy(sim.temperatures, sim.temperatures, sim.potentials, nSteps)
