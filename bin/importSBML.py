#This script will import an SBML model (eg from GeneNetWeaver) and simulate it.
import sys

#install tellurium: https://tellurium.readthedocs.io/en/latest/installation.html
import tellurium as te

print(te.getVersionInfo())

# Load the SBML model from a file
sbml_file = sys.argv[1]
model = te.loadSBMLModel(sbml_file)

# Simulate the model
result = model.simulate(0, 2000, 1000)  # Simulate from time 0 to 100 with 1000 points

# Plot the results
model.plot(result)

#If you want to export file to matlab , that is also possible
#model.exportToMatlab(matlabFilePath)