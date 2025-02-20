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
filename = 'i-C4H10_data.csv'
guessFile = pandas.read_csv(filename)
temperatures = guessFile["T (K)"].tolist()
pressures = guessFile["P (Mpa)"].tolist()
for i in range(len(pressures)):
    pressures[i]*=1E6
structures = guessFile["Structure"].tolist()

'''
compoundData = [0, 0, float(input("Critical Temperature(K): ")),
                float(input("Critical Pressure(MPa): ")),
                float(input("Accentric Factor: ")),
                float(input("PRSK κ1 Value (Optional): ")), 0,
                input("Chemgroup: "),
                float(input("Ionization Potential (eV): ")),
                float(input("Polarizability (Å^3): "))]
H = [float(input("1st Henry's Law Parameter: ")),
     float(input("2nd Henry's Law Parameter: ")),
     float(input("3rd Henry's Law Parameter: ")),
     float(input("4th Henry's Law Parameter: "))]
#TODO: Figure out a better way to do this
pvapConsts = [float(input("1st Empty Hydrate Vapor Pressure Constant: ")),
     float(input("2nd Empty Hydrate Vapor Pressure Constant: ")),
     float(input("3rd Empty Hydrate Vapor Pressure Constant: ")),
     float(input("4th Empty Hydrate Vapor Pressure Constant: "))]
'''

compoundData = [[6, "i-C4H10", float(407.81), float(3.629), float(0.184),
                float(0.0885), 0,
                numpy.array(['{1: 3, 3: 1}'], dtype=object),
                float(10.57),
                float(8.14)]]

H = [190.982, -4913, -34.5102, 0]

pvapConsts = [4.68180,
     -5455.2664,
     2.778907444,
     -0.0089678]

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

def getLangConst(T, P, compoundData, structure):
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
    
    #Water Fugacity Calculation
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
    
    #Get dH from water fugacity
    N_A = 6.022E23
    if structure == "I":
        Vm_hydrate = (11.835+2.217E-5*T+2.242E-6*T**2)**3*(1E-30*N_A/46)-8.006E-9*P/1E6+5.448E-12*(P/1E6)**2
    elif structure == "II":
        Vm_hydrate = (17.13+2.249E-4*T+2.013E-6*T**2-1.009E-9*T**3)**3*(1E-30*N_A/136)-8.006E-9*P/1E6+5.448E-12*(P/1E6)**2
    Psat_hydrate = math.exp(pvapConsts[0]*math.log(T)+pvapConsts[1]/T+pvapConsts[2]+pvapConsts[3]*T)
    
    dH = math.log(f_w/(Psat_hydrate*math.exp(Vm_hydrate*(P-Psat_hydrate)/(R*T))))
    
    if structure == "I":
        nu1 = 0.043478261
        nu2 = 0.130434783
    if structure == "II":
        nu1 = 0.117647059
        nu2 = 0.058823529
    
    frac = [[0],[0]]
    
    def f(frac1):
        if frac1 >=1:
            frac1 = 0.9999999
        frac0 = 0 #1-math.exp((dH-nu2*math.log(1-frac1))/nu1)
        dHexpected = nu1*math.log(1-frac0) + nu2*math.log(1-frac1)
        return abs(dH - dHexpected)
    
    frac[1][0] = scipy.optimize.minimize(f, 0.984, bounds=[(0,0.9999999)]).x
    frac[0][0] = 0 #1-math.exp((dH-nu2*math.log(1-frac[1][0]))/nu1)
    
    Cgg = simFunctions.Lang_GG_Const(T, compoundData, frac, structure)
    Cml = [0,0]
    for i in range(2):
        Cml[i] = frac[i][0]/(Cgg[0][i]*fug_vap*(1-frac[i][0]))
        if Cml [i] < 0:
            Cml [i] = 0
    
    return Cml

def generateParameters(T, lang_consts):
    def model(x, A, B, D):
        return numpy.exp(A+B/x+D/x/x)
    
    initialGuess = [-19.6865, 870.0524, -14208.0553] #Completely arbitrary, based on existing values
    
    A, B, D = scipy.optimize.curve_fit(model, numpy.array(T).flatten(), numpy.array(lang_consts).flatten(), p0=initialGuess)[0]
    return A, B, D

#Logic Begins Here
Cml = [[0 for i in range(len(temperatures))],[0 for i in range(len(temperatures))]]
for i in range(len(temperatures)):
    consts = getLangConst(temperatures[i], pressures[i], compoundData, structures[i])
    Cml[0][i] = consts[0]
    Cml[1][i] = consts[1][0]

A1, B1, D1 = -100, 0, 0 #generateParameters(temperatures, Cml[:][0])
A2, B2, D2 = generateParameters(temperatures, Cml[:][1])

print("Langmuir Constant Parameter Fit:")
print("A1: " + str(A1))
print("B1: " + str(B1))
print("D1: " + str(D1))
print("A2: " + str(A2))
print("B2: " + str(B2))
print("D2: " + str(D2))