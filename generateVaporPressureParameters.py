import simFunctions
import pandas
import numpy
import math
import scipy
from thermo.unifac import UNIFAC, PSRKSG, PSRKIP
import ast

R = 8.31446261815324
errorMargin = 1E-9

fluidProperties = pandas.read_excel('Data.xlsx', sheet_name='Fluid Properties')

#filename = input("Equilibrium Data File name: ")
filename = 'O2 Data.csv'
guessFile = pandas.read_csv(filename)
temperatures = guessFile["T (K)"].tolist()
pressures = guessFile["P (Mpa)"].tolist()
for i in range(len(pressures)):
    pressures[i]*=1E6
structures = guessFile["Structure"].tolist()
smallFrac = guessFile["Small Cage"].tolist()
largeFrac = guessFile["Large Cage"].tolist()

compoundData = [[9, "O2", float(154.581), float(5.043), float(.0222),
                float(0.01512), 0,
                numpy.array(['{119: 1}'], dtype=object),
                float(12.07),
                float(1.562), float(32.00)]]

H = [-286.942, 15450.6, 36.5593, 0.0187662]

LangmuirParameters = [[-121.4476378, 55833.81589, -7364428.344], [-115.7240269, 54879.37555, -7182666.353]]
'''

compoundData = [[10, "H2S", float(373.1), float(9), float(0.1005),
                float(0), 0,
                numpy.array(['{114: 1}'], dtype=object),
                float(10.457),
                float(3.631), float(34.082)]]

H = [-297.158, 16347.7, 40.2024, 0.00257153]

LangmuirParameters = [[-26.33120006, 3822.514824, -45381.84691], [-26.92192792, 7358.818357, -237581.8488]]'''

def Z(compoundData, T, P):
    waterData = numpy.array(fluidProperties.loc[fluidProperties['Compound ID'] == 0])[0]
    localCompoundData = compoundData
    localCompoundData = numpy.column_stack((localCompoundData[0], waterData))
    localCompoundData = localCompoundData.transpose().tolist()
    
    localCompoundData[0][2] = float(localCompoundData[0][2])
    localCompoundData[0][3] = float(localCompoundData[0][3])
    localCompoundData[0][4] = float(localCompoundData[0][4])
    localCompoundData[0][5] = float(localCompoundData[0][5])
    
    Z = simFunctions.PengRobinson(localCompoundData, [0.00001, 0.99999], T, P)[3]
    return Z

def liqPhaseComposition(T, fug_vap, compoundData, P, Psat):
    x = [0, 0]
    Z_inf = Z(compoundData, T, P)
    HenrysLawConst = 101325*math.exp(-1*(H[0] + H[1]/T + H[2]*math.log(T) + H[3]*T))
    x[1] = fug_vap/(HenrysLawConst*math.exp(Z_inf))
    x[0] = 1-sum(x) #Water composition
    return x

def activityCoeff(T, phaseComposition, chemGroups):
    chemGroupList = [{16: 1}]
    for i in range(len(chemGroups)):
        chemGroupList.append(ast.literal_eval(chemGroups[i]))

    GE = UNIFAC.from_subgroups(T, phaseComposition, chemGroupList, interaction_data=PSRKIP, subgroups=PSRKSG).gammas()[0]
    return GE

def freezingPointDepression(T, fug_vap, compoundData, P, chemGroups):
    if T > 273.15:
        Psat = math.exp(4.1539*math.log(T)-5500.9332/T+7.6537-16.1277E-3*T)
    else:
        Psat = math.exp(4.6056*math.log(T)-5501.1243/T+2.9446-T*8.1431E-3)
    
    phaseComposition = liqPhaseComposition(T, fug_vap, compoundData, P, Psat)
    deltadT = R*(273.15)**2/6011*math.log(liqPhaseComposition(T, fug_vap, compoundData, P, Psat)[0]*activityCoeff(T, phaseComposition, chemGroups))
    return deltadT

psatout = numpy.array([0 for i in range(len(temperatures))])

def GetHydratePVap(T, P, fractions, i):
    structure = structures[i]
    
    fug_vap = simFunctions.PengRobinson(compoundData, [1], T, P)[2][0]
    
    if T > 260 and T < 280:
        freezingPoint = 273.15+freezingPointDepression(T, fug_vap, compoundData, P, compoundData[0][7])
        if T > freezingPoint:
            phase = "liquid"
        else:
            phase = "ice"
    elif T <= 260:
        phase = "ice"
    elif T >= 280:
        phase = "liquid"
    
    if phase == "ice":
        Vm_water = 1.912E-5 + T*8.387E-10 + (T**2)*4.016E-12
        Psat_water = math.exp(4.6056*math.log(T)-5501.1243/T+2.9446-T*8.1431E-3)
        f_w = Psat_water*math.exp(Vm_water*(P-Psat_water)/(R*T))
    elif phase == "liquid":
        chemGroups = compoundData[0][7]
        Vm_water = math.exp(-10.921 + 2.5E-4*(T-273.15) - 3.532E-4*(P/1E6-0.101325) + 1.559E-7*((P/1E6-.101325)**2))
        Psat_water = math.exp(4.1539*math.log(T)-5500.9332/T+7.6537-16.1277E-3*T)
        phaseComposition = liqPhaseComposition(T, fug_vap, compoundData, P, Psat_water)
        f_w = phaseComposition[0]*activityCoeff(T, phaseComposition, chemGroups)*Psat_water*math.exp(Vm_water*(P-Psat_water)/(R*T))
        

    cellProperties = simFunctions.getHydrateCellProperties(structure) 
    Deltamu_H_w = 0
    
    Deltamu_H_w += cellProperties[0][10]*math.log(1-fractions[0][i])
    Deltamu_H_w += cellProperties[1][10]*math.log(1-fractions[1][i])
    
    def func(Psat_hydrateGuest):
        return f_w - Psat_hydrateGuest*math.exp(Vm_water*(P-Psat_hydrateGuest)/(R*T))*math.exp(Deltamu_H_w)
    
    Psat = scipy.optimize.fsolve(func, [1000])
    psatout[i] = float(Psat)
    return Psat
    
def generateParameters(T, PVaps):
    def model(x, A, B, D):
        A = abs(A)
        B = -1*abs(B)
        D = -1*abs(D)
        return numpy.exp(A*numpy.log(x)+B/x+2.778907444+D*x)
    
    initialGuess = [4.6446, -5150.369, -0.0087553]
    A, B, D = scipy.optimize.curve_fit(model, numpy.array(T).flatten(), numpy.array(PVaps).flatten(), p0=initialGuess)[0]
    return A, B, D

T = temperatures[i]
P = pressures[i]
fracs= [smallFrac,largeFrac]
pVaps = [0 for i in range(len(temperatures))]
for i in range(len(temperatures)):
    pVaps[i] = float(GetHydratePVap(temperatures[i], pressures[i], fracs, i)[0])
    
pVaps = pVaps
    
A, B, D = generateParameters(temperatures, pVaps)
print("A: " + str(A))
print("B: " + str(B))
print("D: " + str(D))