#This work is licensed under the Creative Commons Attribution-ShareAlike 4.0 International License. To view a copy of this license, visit https://creativecommons.org/licenses/by-sa/4.0/.
#Created by Karsten Kunneman and Amadeu K. Sum at the Colorado School of Mines
#©2025, All Rights Reserved

import simFunctions
import numpy
import time
import math
import csv

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
        logP = numpy.arange(math.log(maxGuessPressure), math.log(minGuessPressure)-(math.log(maxGuessPressure)-math.log(minGuessPressure))/noPoints, -1*(math.log(maxGuessPressure)-math.log(minGuessPressure))/(noPoints-1))
        P = [0 for i in range(len(T))]
        for i in range(len(logP)):
            P[i] = round(math.exp(logP[i]), 2)
    else:
        T = [float(input("Temperature (K): "))]
        P = [float(input("Guess Pressure (MPa): "))*1E6]

exportFileName = input("Export File Name (do not include .csv): ") + ".csv"

numberOfCompounds = int(input("No. of Compounds Excluding Water: "))
components = []
moleFractions = []

IDs, compounds = simFunctions.getComponents()

#Prints and allows user to select from components present in data
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

freshWater = input("Fresh Water? (Y/N): ")
#Prints and allows user to select from salts present in data
salts, inhibitors = simFunctions.getInhibitors()
if freshWater != "Y":
    saltConcs = [0 for i in range(len(salts))]
    inhibitorConcs = [0 for i in range(len(inhibitors))]
    for i in range(len(salts)):
        saltConcs[i] = float(input("Concentration of " + salts[i]+ " (%): "))
    for i in range(len(inhibitors)):
        inhibitorConcs[i] = float(input("Concentration of " + inhibitors[i]+ " (%): "))
    if simFunctions.checkMaxConc(inhibitorConcs) != "":
        raise Exception("Inhibitor(s) " + simFunctions.checkMaxConc(inhibitorConcs) + "Exceed(s) Maximum Concentration")
    if sum(inhibitorConcs)+sum(saltConcs) >= 1:
        raise Exception("Weight Percent of Inhibitors and Salts Exceeds 100%")
        
    betaGas = float(input("betaGas: "))
else:
    saltConcs = [0 for i in range(len(salts))]
    inhibitorConcs = [0 for i in range(len(inhibitors))]
        
print("Calculating...")
startTime = time.time()
eqPressure = numpy.array([0 for i in range(len(T))],dtype=float)
eqStructure = [0 for i in range(len(T))]
eqOccupancy = [0 for i in range(len(T))],[0 for i in range(len(T))]

for i in range(len(T)):
    convergence = simFunctions.equilibriumPressure(T[i], P[i], components, moleFractions, saltConcs, inhibitorConcs)
    eqPressure[i] = convergence[0]/1E6 #In MPa
    eqStructure[i] = convergence[1]
    eqOccupancy[0][i] = convergence[2][0]
    eqOccupancy[1][i] = convergence[2][1]
    print("Temperature " + str(i + 1) + " convergence point reached with a Structure " + convergence[1] + " hydrate.")
    print("Occupancy: " + str(convergence[2]))
    
#Calculate Inhibited Temperatures
betaGas = simFunctions.betaGas(T, eqPressure)
if betaGas == 0:
    betaGas = float(input("betaGas calculation failed, please input new value (leave 0 to skip): "))
TInhibited = [0 for i in range(len(T))]
if freshWater != "Y":
    betaGas = simFunctions.betaGas(T, eqPressure)
    for i in range(len(T)):
        TInhibited[i] = simFunctions.HuLeeSum(T[i], saltConcs, inhibitorConcs, betaGas)
else:
    TInhibited = T
    
endTime = time.time()

print("Time to Complete Calculation: " + str(round(endTime-startTime, 3)) + " seconds")
print("per Data Point: " + str(round((endTime-startTime)/noPoints, 3)) + " seconds")

with open(exportFileName, mode = 'w', newline='') as file:
    data = [['T (K)', 'T Inhibited (K) (If Applicable)', 'P (MPa)', 'Structure', 'Small Cage Occupancy', 'Large Cage Occupancy']]
    for i in range(noPoints):
        data.append([T[i], TInhibited[i], eqPressure[i], eqStructure[i], eqOccupancy[0][i], eqOccupancy[1][i]])
    writer = csv.writer(file)
    writer.writerows(data)