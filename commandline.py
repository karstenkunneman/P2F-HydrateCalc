import simFunctions
import numpy
import time
import matplotlib.pyplot as plt
import math

usecsv = input("Import .csv file? (Y/N): ")
if usecsv == "Y":
    filename = input("File name: ")
    guessFile = numpy.genfromtxt(filename, delimiter=',', skip_header=1)
    T = guessFile[:,0]
    P = guessFile[:,1]*1E6
    noPoints = len(T)
else:
    noPoints = int(input("Number of Data Points: "))
    if noPoints > 1:
        minTemp = float(input("Lower Temperature (K): ")) #Temp in K
        maxTemp = float(input("Upper Temperature (K): "))
        minGuessPressure = float(input("Minimum Guess Pressure (MPa): "))*1E6 #Pressure in Pa
        maxGuessPressure = float(input("Maximum Guess Pressure (MPa): "))*1E6 #Pressure in Pa
        T = numpy.arange(maxTemp, minTemp-(maxTemp-minTemp)/noPoints, -1*(maxTemp-minTemp)/(noPoints-1))
        P = numpy.arange(maxGuessPressure, minGuessPressure-(maxGuessPressure-minGuessPressure)/noPoints, -1*(maxGuessPressure-minGuessPressure)/(noPoints-1))
    else:
        T = [float(input("Temperature (K): "))]
        P = [float(input("Guess Pressure (MPa): "))*1E6]
    
numberOfCompounds = int(input("No. of Compounds Excluding Water: "))
components = []
moleFractions = []

IDs, compounds = simFunctions.getComponents()

for i in range(math.ceil(len(IDs)/2)):
    try:
        if len(str(IDs[2*i]) + ". " + compounds[2*i]) < 7:
            print(" " + str(IDs[2*i]) + ". " + compounds[2*i] + "\t" + "\t" + str(IDs[2*i+1]) + ". " + compounds[2*i+1])
        else:
            print(" " + str(IDs[2*i]) + ". " + compounds[2*i] + "\t" + str(IDs[2*i+1]) + ". " + compounds[2*i+1])
    except:
        print(" " + str(IDs[2*i]) + ". " + compounds[2*i])
    
for i in range(numberOfCompounds):
    components += [int(input("Compound ID " + str(i + 1) + " : "))]
    if numberOfCompounds > 1:
        if i != numberOfCompounds-1:
            moleFractions += [float(input("Mole Fraction of Compound " + str(i + 1) + " : "))]
        else:
            moleFractions.append(1-sum(moleFractions))
    else:
        moleFractions = [1]

print("Calculating...")
startTime = time.time()
eqPressure = numpy.array([0 for i in range(len(T))],dtype=float)
eqStructure = [0 for i in range(len(T))]

for i in range(len(T)):
    convergence = simFunctions.equilibriumPressure(T[i], P[i], components, moleFractions)
    eqPressure[i] = convergence[0]/1E6 #In MPa
    eqStructure[i] = convergence[1]
    print("Temperature " + str(i + 1) + " convergence point reached with a Structure " + convergence[1] + " hydrate.")
    print("Occupancy: " + str(convergence[2]))

plt.plot(T, eqPressure, '-ok')
plt.yscale("log")
plt.title("Equilibrium Predictions")
plt.xlabel("Temperature (K)")
plt.ylabel("Pressure (MPa)")
plt.show

endTime = time.time()

print("Time to Complete Calculation: " + str(round(endTime-startTime, 3)) + " seconds")
print("per Data Point: " + str(round((endTime-startTime)/noPoints, 3)) + " seconds")