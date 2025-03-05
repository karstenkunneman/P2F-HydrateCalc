import simFunctions
import pandas

filename = 'CH4.csv'
guessFile = pandas.read_csv(filename)

T = guessFile["T (K)"].tolist()
eqPressure = guessFile["P (Mpa)"].tolist()
eqStructure = "I"

beta = simFunctions.betaGas(T, eqPressure, eqStructure)

print(beta)