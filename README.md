# Barataria Basin Sediment Budget Model

Model that implements the Monte Carlo sampling and solution of equation 1 in Edmonds et al.

Barataria_SedimentBudget.m starts a series of 5 Monte Carlo simulations for the different scenarios described in Edmonds et al.

Barataria_SedimentBudget.m calls function BaratariaArea_final.m to calculate how the Barataria's area changes for each set of parameters selected via the Monte Carlo. 

CalculateOverbankDeposition.m ingests the Delft3D output files and calculates the scaled-up yearly mass deposition rate in Barataria. 
