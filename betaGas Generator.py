import pandas
import numpy
import math
import matplotlib.pyplot as plt

R = 8.31446261815324

fluidProperties = pandas.read_excel('Data.xlsx', sheet_name='Fluid Properties')

filename = input("Equilibrium Data File name: ")
guessFile = pandas.read_csv(filename)
temperatures = guessFile["T (K)"].tolist()
pressures = guessFile["P (Mpa)"].tolist()
for i in range(len(pressures)):
    pressures[i]*=1E6

inverseTemp = [0 for i in range(len(temperatures))]
for i in range(len(temperatures)):
    inverseTemp[i] = 1/temperatures[i]
    
lnPressure = [0 for i in range(len(temperatures))]
for i in range(len(temperatures)):
    lnPressure[i] = math.log(pressures[i])
    
for i in range(len(temperatures)):
    if temperatures[i] <= 273.15:
        del lnPressure[i]
        
inverseTemp = [x for x in inverseTemp if x <= 1/273.15]
    
slope = numpy.polyfit(inverseTemp, lnPressure, 1)[0]

betaGas = -5.75/slope

print(betaGas)

plt.plot(inverseTemp, lnPressure, '-ok')
plt.yscale("log")
plt.title("Equilibrium Predictions")
plt.xlabel("Temperature (K)")
plt.ylabel("Pressure (MPa)")
plt.show