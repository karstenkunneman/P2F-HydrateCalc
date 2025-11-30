#This work is licensed under the Creative Commons Attribution-ShareAlike 4.0 International License. To view a copy of this license, visit https://creativecommons.org/licenses/by-sa/4.0/.
#Created by Karsten Kunneman and Amadeu K. Sum at the Colorado School of Mines
#©2025, All Rights Reserved

import numpy
import scipy
import csv
import p2f_HydrateCalcLib.core as core
import p2f_HydrateCalcLib.model as model

#REMOVE EVENTUALLY
import warnings
warnings.filterwarnings("ignore")

#Constants
errorMargin = 1E-4 #tested to 1 Pa accuracy difference
R = 8.31446261815324 #m^3 Pa/mol K
boltzmannConstant = 1.380649E-23 #m^2 kg/s^2 K
N_A = 6.022E23

def getComponents():
    IDs = core.fluidData[:, 0].tolist()
    del IDs[0]
    compounds = core.fluidData[:, 1].tolist()
    del compounds[0]
    
    return IDs, compounds

def getInhibitors():
    salts = core.saltData[:,0]
    inhibitors = core.inhibitorData[:,0]
    return salts, inhibitors

#For a given hydrate structure, return its radii and coordination numbers
def getHydrateCellProperties(structure):
    cellProperties = numpy.array(model.hydrateCellProperties.loc[(model.hydrateCellProperties['Structure'] == structure)])
    return cellProperties

#Ensure that the concentration of no organic inhibitor exceeds its individual maximum concentration
def checkMaxConc(inhibitorConcs):
    exceededInhibitors = ""
    for i in range(len(inhibitorConcs)):
        if inhibitorConcs[i] > core.inhibitorData[i][5]:
            exceededInhibitors += str(core.inhibitorData[i][0]) + " "
    return exceededInhibitors

#For a given salt or inhibitor, determine its minimum concentration to achieve a given temperature inhibition from a given pure water equilibrium temperature
def getConcentration(T, TDesired, inhibitor, salt, betaGas, freezingPoint):
    saltList, inhibitorList = getInhibitors()
    
    inhibitorConcs = numpy.zeros(len(inhibitorList))
    saltConcs = numpy.zeros(len(saltList))

    def f(conc, inhibitor):
        if inhibitor != "salt":
            for i in range(len(inhibitorList)):
                if i == inhibitor:
                    inhibitorConcs[i] = conc
                else:
                    inhibitorConcs[i] = 0
        else:
            for i in range(len(saltList)):
                if i == salt:
                    saltConcs[i] = conc
                else:
                    saltConcs[i] = 0
        Tinhibited = core.HuLeeSum(T, saltConcs, inhibitorConcs, betaGas, freezingPoint)
        return TDesired-Tinhibited
        
    conc = scipy.optimize.fsolve(f,1,args=inhibitor)[0]

    return conc

def tempConversion(tempUnit, T, toK):
    if toK == False:
        if tempUnit == "K":
            T = T
        elif tempUnit == "°C":
            T = T - 273.15
        elif tempUnit == "°F":
            T = (T-273.15)*1.8 + 32
    else:
        if tempUnit == "K":
            T = T
        elif tempUnit == "°C":
            T = T + 273.15
        elif tempUnit == "°F":
            T = (T-32)/1.8 + 273.15
    return T
    
def pressureConversion(pressureUnit, P, toMPa):
    if toMPa == False:
        if pressureUnit == "MPa":
            P = P
        elif pressureUnit == "psia":
            P *= 145.038
        elif pressureUnit == "bar":
            P *= 10
    else: 
        if pressureUnit == "MPa":
            P = P
        elif pressureUnit == "psia":
            P /= 145.038
        elif pressureUnit == "bar":
            P /= 10
    return P

def guestComp(thetaSmall, thetaLarge, structure):
    normalizedComps = numpy.zeros(len(thetaSmall))

    for i in range(len(normalizedComps)):
        if structure == "I":
                normalizedComps[i] = 2*thetaSmall[i]+6*thetaLarge[i]
        else:
                normalizedComps[i] = 16*thetaSmall[i]+8*thetaLarge[i]

    fracs = numpy.zeros(len(thetaSmall))
    sumComps = sum(normalizedComps)
    for i in range(len(normalizedComps)):
        if structure == "I":
                fracs[i] = round((2*thetaSmall[i]+6*thetaLarge[i])/sumComps,4)
        else:
                fracs[i] = round((16*thetaSmall[i]+8*thetaLarge[i])/sumComps,4)

    return fracs

def generateOutput(componentNames, componentIDs, moleFractions, salts, saltConcs, inhibitors, 
                   inhibitorConcs, T, TInhibited, P, convergence, IDs, tempUnit, pressureUnit):
    numpy.set_printoptions(suppress=True)
    with open("Output Template.csv", newline='', encoding='utf-8-sig') as csvfile:
        reader = list(csv.reader(csvfile))

        if not isinstance(moleFractions[0], list):
            moleFractions = [moleFractions*100]
            componentIDs = [componentIDs]

        #Add salts and inhibitors
        inhibitorLineNo = 0
        for i in range(len(salts)):
            if saltConcs[i] != 0:
                reader.insert(3+inhibitorLineNo, [salts[i],saltConcs[i], None, None, None])
                inhibitorLineNo += 1
        for i in range(len(inhibitors)):
            if inhibitorConcs[i] != 0:
                reader.insert(3+inhibitorLineNo, [inhibitors[i],inhibitorConcs[i], None, None, None])
                inhibitorLineNo += 1
               
        try:
            if len(set(tuple(row) for row in moleFractions)) == 1 and len(set(tuple(row) for row in componentIDs)) == 1:
                j = 0
                for i in range(len(IDs)):
                    if int(i+1) in componentIDs[0]:
                        if j < inhibitorLineNo:
                            reader[3+j][3] = componentNames[i]
                            reader[3+j][4] = round(moleFractions[0][j]*100,2)
                            j += 1
                        else:
                            reader.insert(3+j, [None, None, None, componentNames[i],moleFractions[0][j]*100])
                            j += 1

                reader = [row[:11] for row in reader]
        except:
            j = 0
            for i in range(len(IDs)):
                if int(i+1) in componentIDs:
                    if j < inhibitorLineNo:
                        reader[3+j][3] = componentNames[i]
                        reader[3+j][4] = moleFractions[0][j]
                        j += 1
                    else:
                        reader.insert(3+j, [None, None, None, componentNames[i],moleFractions[j]*100])
                        j += 1

            reader = [row[:10] for row in reader]

        headerRow = len(reader)
        reader[headerRow-1][0] = "T (" + tempUnit + ")"
        reader[headerRow-1][1] = "Inhib. T (" + tempUnit + ")"
        reader[headerRow-1][2] = "P (" + pressureUnit + ")"

        for i in range(len(convergence)):
            try:
                insertList = [T[i], TInhibited[i], P[i], convergence[i][2], str([round(value, 4) for value in convergence[i][3][0]]), str([round(value, 4) for value in convergence[i][3][1]]), round(convergence[i][4],2), round(convergence[i][5][0],1), round(convergence[i][5][1],1), convergence[i][6], str(guestComp(convergence[i][3][0].tolist(), convergence[i][3][1].tolist(), convergence[i][2]).tolist())]
            except:
                insertList = [T[i], None, P[i], convergence[i][2], str([round(value, 4) for value in convergence[i][3][0]]), str([round(value, 4) for value in convergence[i][3][1]]), round(convergence[i][4],2), round(convergence[i][5][0],1), round(convergence[i][5][1],1), convergence[i][6], str(guestComp(convergence[i][3][0].tolist(), convergence[i][3][1].tolist(), convergence[i][2]).tolist())]
            try:
                if len(set(tuple(row) for row in moleFractions)) != 1 or len(set(tuple(row) for row in componentIDs)) != 1:
                    for j in range(len(IDs)):
                        if int(j+1) in componentIDs[i]:
                            index = componentIDs[i].index(int(j+1))
                            insertList.append(moleFractions[i][index]*100)
                        else:
                            insertList.append('')
                    del(reader[2][3:])
            except:
                for j in range(len(IDs)):
                    if int(j+1) in componentIDs:
                        index = componentIDs.index(int(j+1))
                        insertList.append(moleFractions[index]*100)
                    else:
                        insertList.append('')
                del(reader[2][3:])
            reader.insert(5+len(salts)+len(inhibitors)+i, insertList)
        
    return reader

def massToMolFrac(compounds, weightFracs):
    compoundData = numpy.array(core.fluidData[core.fluidData[:, 1] == compounds[0]])
    for i in range(len(compounds) - 1):
        compoundData = numpy.vstack([compoundData, numpy.array(core.fluidData[core.fluidData[:, 1] == compounds[i+1]])])

    molWeights = []
    for i in range(len(compounds)):
        molWeights.append(compoundData[i][10])

    mols = []
    for i in range(len(compounds)):
        mols.append(weightFracs[i]/molWeights[i])

    moleFracs = []
    for i in range(len(compounds)):
        try:
            moleFracs.append(mols[i]/sum(mols)*sum(weightFracs))
        except:
            moleFracs.append(0)
        
    return moleFracs