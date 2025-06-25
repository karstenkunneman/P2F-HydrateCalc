#This work is licensed under the Creative Commons Attribution-ShareAlike 4.0 International License. To view a copy of this license, visit https://creativecommons.org/licenses/by-sa/4.0/.
#Created by Karsten Kunneman and Amadeu K. Sum at the Colorado School of Mines
#©2025, All Rights Reserved

import math 
import numpy
import pandas
import scipy
from thermo.unifac import UNIFAC, PSRKSG, PSRKIP
import ast

#REMOVE EVENTUALLY
import warnings
warnings.filterwarnings("ignore")

#Extract data from excel files
fluidProperties = pandas.read_excel('Data.xlsx', sheet_name='Fluid Properties')
guessConstants = pandas.read_excel('Data.xlsx', sheet_name='Guess Constants')
hydrateCellProperties = pandas.read_excel('Data.xlsx', sheet_name='Hydrate Cell Properties')
vaporPressureConstants = pandas.read_excel('Data.xlsx', sheet_name='Vapor Pressure Constants')
mixConstants = pandas.read_excel('Data.xlsx', sheet_name='Binary Interaction Parameters')
saltData = pandas.read_excel('Data.xlsx', sheet_name='Salts Data').to_numpy()
inhibitorData = pandas.read_excel('Data.xlsx', sheet_name='Inhibitor Data').to_numpy()
langParameters = pandas.read_excel('Data.xlsx', sheet_name='Langmuir Parameters')
henrysLawConstants = pandas.read_excel('Data.xlsx', sheet_name='Henrys Law Parameters')

#Constants
errorMargin = 1E-4 #tested to 1 Pa accuracy difference
R = 8.31446261815324 #m^3 Pa/mol K
boltzmannConstant = 1.380649E-23 #m^2 kg/s^2 K

def getComponents():
    IDs = fluidProperties["Compound ID"].tolist()
    del IDs[0]
    compounds = fluidProperties["Formula"].tolist()
    del compounds[0]
    
    return IDs, compounds

def getInhibitors():
    salts = saltData[:,0]
    inhibitors = inhibitorData[:,0]
    return salts, inhibitors

#PRSV Equation of State
def PengRobinson(compoundData, moleFractions, T, P):
    noComponents = len(moleFractions)
    #Calculate "b" value for each pure substance and add to weighted and unweighted lists
    bmix = 0
    b = [0]*noComponents
    for i in range(noComponents):    
        Tc = compoundData[i][2]
        Pc = compoundData[i][3]*1E6
        b[i] = 0.07780*(R*Tc/Pc)
        bmix += b[i]*moleFractions[i]
        
    #Calculate the "a" values for all pure compounds
    a = [0]*noComponents
    for i in range(noComponents):
        Tc = compoundData[i][2]
        Pc = compoundData[i][3]*1E6
        omega = compoundData[i][4]
        kappa1 = compoundData[i][5]
        kappa = 0.378893 + 1.4897153*omega - 0.17131848*omega**2 + 0.0196554*omega**3 + kappa1*(1+math.sqrt(T/Tc))*(0.7-(T/Tc))
        #kappa = 0.37464 + 1.54226*omega - 0.26992*(omega**2)
        alpha = (1 + kappa*(1-math.sqrt(T/Tc)))**2
        a[i] = 0.45724*((R**2)*(Tc**2)/Pc)*alpha
        
    #Obtain the interaction parameters for all combinations of compounds
    if noComponents > 1:
        interactionParameters = numpy.zeros((noComponents,noComponents))
        for i in range(noComponents):
            for j in range(noComponents):
                interactionParameters[i][j] = mixConstants.loc[(mixConstants['Compound 1'] == compoundData[i][1]), (compoundData[j][1])].reset_index(drop=True)[0]
    else:
        interactionParameters = [[0]]
        
    #Finally, determine the partial "a" values based on mole ratios
    amix = 0
    xia = numpy.zeros(noComponents)
    for i in range(noComponents):
        for j in range(noComponents):
            amix += (math.sqrt(a[i]*a[j]))*(1-interactionParameters[i][j])*moleFractions[i]*moleFractions[j]
            xia[i] += (math.sqrt(a[i]*a[j]))*(1-interactionParameters[i][j])*moleFractions[j]
            
    A = amix*P/((R*T)**2)
    B = (bmix*P)/(R*T)
    
    Z = numpy.roots([1, -1+B, A-3*(B**2)-2*B, -1*A*B+B**2+B**3])
    
    #Analyze roots
    noReal = 0
    realZ = []
    for i in range(len(Z)):
        if numpy.iscomplex(Z[i]) == False:
            noReal += 1
            realZ.append(Z[i])
    
    if noReal == 3:
        ZVmix = max(realZ).real
        ZLmix = min(realZ).real
    elif noReal == 1:
        ZVmix = realZ[0].real
        ZLmix = realZ[0].real
    
    #Calculate molar volumes for each compound
    VmVap = (ZVmix*R*T/P)
    VmLiq = (ZLmix*R*T/P)
    
    fugVap = numpy.zeros(noComponents)
    
    #Calculate fugacities for each of the compounds
    for i in range(noComponents):
        fugVap[i] = P*moleFractions[i]*math.exp((b[i]/bmix)*(ZVmix-1)-math.log(ZVmix-B)-(A/(2*math.sqrt(2)*B))*(2*xia[i]/amix-b[i]/bmix)*math.log((ZVmix+(1+math.sqrt(2))*B)/(ZVmix+(1-math.sqrt(2))*B)))
        #fugLiq[i] = P*moleFractions[i]*math.exp((b[i]/bmix)*(ZLmix-1)-math.log(ZLmix-B)-(A/(2*math.sqrt(2)*B))*(2*xia[i]/amix-b[i]/bmix)*math.log((ZLmix+(1+math.sqrt(2))*B)/(ZLmix+(1-math.sqrt(2))*B)))
       
    return VmVap, VmLiq, fugVap, ZVmix

#For a given hydrate structure, return its radii and coordination numbers
def getHydrateCellProperties(structure):
    cellProperties = numpy.array(hydrateCellProperties.loc[(hydrateCellProperties['Structure'] == structure)])
    return cellProperties

#Calculates Langmuir Constant Cml
def Lang_Const(T, cellRadii, a, RCell, z, epsilon, sigma, structure, shell, compound):    
    condition = (langParameters['Structure'] == structure) & \
        (langParameters['Shell'] == shell) & \
        (langParameters['Guest'] == compound)
    
    Ac = langParameters.loc[condition, "Ac"].values[0]
    Bc = langParameters.loc[condition, "Bc"].values[0]
    Dc = langParameters.loc[condition, "Dc"].values[0]
    
    Cml = math.exp(Ac + Bc/T + Dc/T/T)
    
    return Cml

#Calculates Langmuir Guest-Guest Parameter Cgg
def Lang_GG_Const(T, compoundData, fracs, structure):
    try:
        I = compoundData[:,8]#eV
        alpha = compoundData[:,9]#angstrom^3
    except:
        I = [compoundData[0][8]]
        alpha = [compoundData[0][9]]
        
    noComponents = len(I)    
        
    Cgg = [[1,1] for i in range(noComponents)]
    
    C6 = numpy.zeros((noComponents,noComponents))
    C8 = numpy.zeros((noComponents,noComponents))
    C10 = numpy.zeros((noComponents,noComponents))
            
    for i in range(noComponents):
        for j in range(noComponents):
            C6[i][j] = (3/2)*alpha[i]*alpha[j]*I[i]*I[j]/(I[i]+I[j])*23.05/.001987
            C8[i][j] = 496778.3824*alpha[i]*alpha[j]*(I[i]/(2*I[i]+I[j])+I[j]/(2*I[j]+I[i]))
            C10[i][j] = 13260598.42*alpha[i]*alpha[j]/(I[i]+I[j])
    
    wrgg = [[[[0,0],[0,0]] for i in range(noComponents)] for i in range(noComponents)]
    
    #Unsure exactly what these are relative to r
    if structure == "I":
        a_lc = 12.03
        for i in range(noComponents):
            for j in range(noComponents):
                r0 = a_lc*0.86602
                wrgg[i][j][0][0] = -1*C6[i][j]*12.25367/(r0**6)-C8[i][j]*10.3552/(r0**8)-C10[i][j]*9.5644/(r0**10)
                r1 = a_lc*0.55901
                wrgg[i][j][0][1] = -1*C6[i][j]*13.41525/(r1**6)-C8[i][j]*12.38994/(r1**8)-C10[i][j]*12.12665/(r1**10)
                
                r0 = a_lc*0.55901
                wrgg[i][j][1][0] = -1*C6[i][j]*4.47175/(r0**6)-C8[i][j]*4.12998/(r0**8)-C10[i][j]*4.02482/(r0**10)
                r1 = a_lc*0.5
                wrgg[i][j][1][1] = -1*C6[i][j]*5.14048/(r1**6)-C8[i][j]*3.74916/(r1**8)-C10[i][j]*3.09581/(r1**10)
    elif structure == "II":
        a_lc = 17.31
        for i in range(noComponents):
            for j in range(noComponents):
                r0 = a_lc*0.35355
                wrgg[i][j][0][0] = -1*C6[i][j]*6.92768/(r0**6)-C8[i][j]*6.23392/(r0**8)-C10[i][j]*6.06724/(r0**10)
                r1 = a_lc*0.41457
                wrgg[i][j][0][1] = -1*C6[i][j]*6.91143/(r1**6)-C8[i][j]*6.28271/(r1**8)-C10[i][j]*6.10214/(r1**10)
                
                r0 = a_lc*0.41457
                wrgg[i][j][1][0] = -1*C6[i][j]*13.82287/(r0**6)-C8[i][j]*12.56542/(r0**8)-C10[i][j]*12.20428/(r0**10)
                r1 = a_lc*0.43301
                wrgg[i][j][1][1] = -1*C6[i][j]*5.11677/(r1**6)-C8[i][j]*4.33181/(r1**8)-C10[i][j]*4.11102/(r1**10)
    
    #Obtain the interaction energy of guestz`
    for i in range(noComponents):
        for j in range(noComponents):
            for k in range(2):
                Cgg[i][k] *= math.exp(-1*(wrgg[i][j][k][0]*fracs[0][j])/T)*math.exp(-1*(wrgg[i][j][k][1]*fracs[1][j])/T)
    
    return Cgg

#Calculates the fractional occupancy of small and large shells by each component
def frac(T, kiharaParameters, vaporFugacities, compoundData, structure, compounds):
    noComponents = len(vaporFugacities)
    cellProperties = getHydrateCellProperties(structure) 
    guessFractions = numpy.zeros((2, noComponents))
    oldGuessFractions = [[1 for i in range(noComponents)],[1 for i in range(noComponents)]]
    Cgg = [[1.5, 1.5] for i in range(noComponents)]
    langConsts = [[None, None] for i in range(noComponents)]
    RCell = [0,0]
    fracDiff = [[1 for i in range(noComponents)],[1 for i in range(noComponents)]]
    while abs(sum(fracDiff[0])/noComponents+sum(fracDiff[1])/noComponents)/2 >= errorMargin: #average fractional occupancy difference of all shells
        fracDiff = numpy.zeros((2, noComponents))
        for i in range(2):
            denominator = 0
            #Hydrate Cell Properties
            cellRadii = [0, 0, 0]
            for j in range(3):
                cellRadii[j] = cellProperties[i][2+j]*1E-10
            z = [0, 0, 0]
            for j in range(3):
                z[j] = cellProperties[i][5+j]
            RCell[i] = cellProperties[i][8]*1E-10
                   
            for j in range(noComponents):
                langConsts[j][i] = Lang_Const(T, cellRadii, None, RCell[i], z, None, None, structure, i, compounds[j])*math.sqrt(1-fracDiff[i][j])
                denominator += Cgg[j][i]*langConsts[j][i]*vaporFugacities[j]
                    
            for j in range(noComponents):
                guessFractions[i][j] = Cgg[j][i]*langConsts[j][i]*vaporFugacities[j]/(1 + denominator)
                fracDiff[i][j] = abs(guessFractions[i][j]/oldGuessFractions[i][j]-1)/noComponents
                oldGuessFractions[i][j]=guessFractions[i][j]
            
        Cgg = Lang_GG_Const(T, compoundData, guessFractions, structure)

    for i in range(2):
        for j in range(noComponents):
            if guessFractions[i][j] < 1E-10:
                guessFractions[i][j] = 0
        
    return guessFractions

#Equation 2
def deltaHydratePotential(T, kiharaParameters, structure, vaporFugacities, compoundData, compounds):
    cellProperties = getHydrateCellProperties(structure) 
    fractions = 0
    Deltamu_H_w = 0

    fractions = frac(T, kiharaParameters, vaporFugacities, compoundData, structure, compounds)
    
    Deltamu_H_w += cellProperties[0][10]*math.log(1-sum(fractions[0]))
    Deltamu_H_w += cellProperties[1][10]*math.log(1-sum(fractions[1]))
    return Deltamu_H_w, fractions

#Equation 18
def henrysLawConst(compound, T):
    constants = numpy.array(henrysLawConstants.loc[henrysLawConstants['Compound ID'] == compound])[0]
    H_i = 101325*math.exp(-1*(constants[2]/1.987 + constants[3]/T/1.987 + constants[4]*math.log(T)/1.987 + constants[5]*T/1.987))
    return H_i

#Infinite Dilution Compressibility Factor
def Z(compoundData, T, P):
    noComponents = len(compoundData)
    Z = numpy.zeros(noComponents)
    waterData = numpy.array([0, "H2O", 647.3, 22.09, 0.3438, -0.06635, None, {16: 1}, None, None, 18.02])
    for i in range(noComponents):
        localCompoundData = compoundData[i]
        localCompoundData = numpy.column_stack((localCompoundData, waterData))
        localCompoundData = localCompoundData.transpose()
        Z[i] = PengRobinson(localCompoundData, [0.00001, 0.99999], 273.15, P)[3]
    return Z

#Equation 17
def liqPhaseComposition(compounds, T, fug_vap, compoundData, P, Psat):
    noComponents = len(compounds)
    x = numpy.zeros(noComponents+1)
    for i in range(noComponents):
        Z_inf = Z(compoundData, T, P)
        x[i+1] = fug_vap[i]/(henrysLawConst(compounds[i], T)*math.exp(Z_inf[i]))
    x[0] = 1-sum(x) #Water composition
    return x

#Equation 19
def freezingPointDepression(compounds, T, fug_vap, compoundData, P, chemGroups, saltConcs, inhibitorConcs):
    if T > 273.15:
        Psat = math.exp(4.1539*math.log(T)-5500.9332/T+7.6537-16.1277E-3*T)
    else:
        Psat = math.exp(4.6056*math.log(T)-5501.1243/T+2.9446-T*8.1431E-3)
    
    phaseComposition = liqPhaseComposition(compounds, T, fug_vap, compoundData, P, Psat)
    deltadT = R*(273.15)**2/6011*math.log(liqPhaseComposition(compounds, T, fug_vap, compoundData, P, Psat)[0]*activityCoeff(T, phaseComposition, chemGroups))
    return deltadT

def activityCoeff(T, phaseComposition, chemGroups):
    chemGroupList = [{16: 1}]
    for i in range(len(chemGroups)):
        chemGroupList.append(ast.literal_eval(chemGroups[i]))

    GE = UNIFAC.from_subgroups(T, phaseComposition, chemGroupList, interaction_data=PSRKIP, subgroups=PSRKSG).gammas()[0]
    return GE
    
#Equations 15 and 16 
def waterFugacity(T, P, phase, fug_vap, compounds, compoundData):
    if phase == "ice":
        Vm_water = 1.912E-5 + T*8.387E-10 + (T**2)*4.016E-12
        Psat_water = math.exp(4.6056*math.log(T)-5501.1243/T+2.9446-T*8.1431E-3)
        f_w = Psat_water*math.exp(Vm_water*(P-Psat_water)/(R*T))
    elif phase == "liquid":
        chemGroups = compoundData[:,7]
        Vm_water = math.exp(-10.921 + 2.5E-4*(T-273.15) - 3.532E-4*(P/1E6-0.101325) + 1.559E-7*((P/1E6-.101325)**2))
        Psat_water = math.exp(4.1539*math.log(T)-5500.9332/T+7.6537-16.1277E-3*T)
        phaseComposition = liqPhaseComposition(compounds, T, fug_vap, compoundData, P, Psat_water)
        f_w = phaseComposition[0]*activityCoeff(T, phaseComposition, chemGroups)*Psat_water*math.exp(Vm_water*(P-Psat_water)/(R*T))
    return f_w

#Equation A2
def hydrateFugacity(T, P, PvapConsts, structure, fug_vap, compounds, kiharaParameters, compoundData):
    noComponents = len(compounds)
    cellProperties = getHydrateCellProperties(structure) 
    N_A = 6.022E23
    if structure == "I":
        Vm_water = (11.835+2.217E-5*T+2.242E-6*T**2)**3*(1E-30*N_A/46)-8.006E-9*P/1E6+5.448E-12*(P/1E6)**2
    elif structure == "II":
        Vm_water = (17.13+2.249E-4*T+2.013E-6*T**2-1.009E-9*T**3)**3*(1E-30*N_A/136)-8.006E-9*P/1E6+5.448E-12*(P/1E6)**2
    
    dH = deltaHydratePotential(T, kiharaParameters, structure, fug_vap, compoundData, compounds)
    frac = dH[1]
    
    A = 0
    B = 0
    D = 0
    
    denominator = 0
    for i in range(noComponents):
        denominator += frac[0][i]*cellProperties[0][10]
        denominator += frac[1][i]*cellProperties[1][10]
            
    Z = numpy.zeros(noComponents)
    for i in range(noComponents):
        for j in range(2):
            Z[i] += (frac[j][i]*cellProperties[j][10])/denominator
            
    for i in range(noComponents):
        A += PvapConsts[i,3]*Z[i]
        B += PvapConsts[i,4]*Z[i]
        D += PvapConsts[i,6]*Z[i]
    
    Psat_water = math.exp(A*math.log(T)+B/T+2.7789+D*T)
    
    f_h = Psat_water*math.exp(Vm_water*(P-Psat_water)/(R*T))*math.exp(dH[0])
    return f_h,frac

def HuLeeSum(T, saltConcs, inhibitorConcs, betaGas):
    noSalts = len(saltConcs)
    noInhibitors = len(inhibitorConcs)
    if T >= 273.15:
        if sum(saltConcs) != 0:
            catMols = numpy.zeros(noSalts)
            anMols = numpy.zeros(noSalts)
            catMolFracs = numpy.zeros(noSalts)
            anMolFracs = numpy.zeros(noSalts)
            for i in range(noSalts):
                catMols[i] = saltConcs[i]/saltData[i][1]*saltData[i][4]
                anMols[i] = saltConcs[i]/saltData[i][1]*saltData[i][5]
            totalSaltMols = sum(catMols)+sum(anMols)
            waterMols = (100-sum(saltConcs))/18.015
            for i in range(noSalts):
                catMolFracs[i] = saltData[i][2]*catMols[i]/(waterMols+totalSaltMols)
                anMolFracs[i] = saltData[i][3]*anMols[i]/(waterMols+totalSaltMols)
                
            inhibitorMols = numpy.zeros(noInhibitors)
            inhibitorMolFracs = numpy.zeros(noInhibitors)
            for i in range(noInhibitors):
                inhibitorMols[i] = inhibitorConcs[i]/(1-0.01*inhibitorConcs[i])/inhibitorData[i][1]
            totalMols = totalSaltMols + sum(inhibitorMols)
            for i in range(noInhibitors):
                inhibitorMolFracs[i] = inhibitorMols[i]/totalMols

            totalEffSaltMols = sum(catMolFracs) + sum(anMolFracs)
        else:
            waterMols = (100-sum(saltConcs))/18.015
            totalSaltMols = 0
            totalEffSaltMols = 0

        inhibitorMols = numpy.zeros(noInhibitors)
        inhibitorMolFracs = numpy.zeros(noInhibitors)
        for i in range(noInhibitors):
            inhibitorMols[i] = inhibitorConcs[i]/(1-0.01*inhibitorConcs[i])/inhibitorData[i][1]
        totalMols = totalSaltMols + sum(inhibitorMols)
        for i in range(noInhibitors):
            inhibitorMolFracs[i] = inhibitorMols[i]/(totalMols+waterMols)

        lnawSalts = -1.06152*totalEffSaltMols+3.25726*totalEffSaltMols*totalEffSaltMols-37.2263*totalEffSaltMols*totalEffSaltMols*totalEffSaltMols
        lnawInhibitors = inhibitorMols = numpy.zeros(noInhibitors)
        for i in range(noInhibitors):
            lnawInhibitors[i] = inhibitorData[i][2]*inhibitorMolFracs[i]+inhibitorData[i][3]*inhibitorMolFracs[i]**2

        Tinhibited = T/(1-betaGas*(lnawSalts+sum(lnawInhibitors))*T)
    else:
        Tinhibited = None
    return Tinhibited

def checkMaxConc(inhibitorConcs):
    exceededInhibitors = ""
    for i in range(len(inhibitorConcs)):
        if inhibitorConcs[i] > inhibitorData[i][5]:
            exceededInhibitors += str(inhibitorData[i][0]) + " "
    return exceededInhibitors

def getConcentration(T, TDesired, inhibitor, salt, betaGas, noInhibitors, noSalts):
    inhibitorConcs = numpy.zeros(noInhibitors)
    saltConcs = numpy.zeros(noSalts)
    
    def f(conc, inhibitor):
        if inhibitor != "salt":
            for i in range(noInhibitors):
                if i == inhibitor:
                    inhibitorConcs[i] = conc[0]
                else:
                    inhibitorConcs[i] = 0
        else:
            for i in range(noSalts):
                if i == salt:
                    saltConcs[i] = conc[0]
                else:
                    saltConcs[i] = 0
        Tinhibited = HuLeeSum(T, saltConcs, inhibitorConcs, betaGas)
        return TDesired-Tinhibited
        
    conc = scipy.optimize.fsolve(f,1,args=inhibitor)[0]

    return conc

def betaGas(temperatures, pressures):
    for i in range(len(pressures)):
        pressures[i]

    inverseTemp = []
    lnPressure = []
    for i in range(len(temperatures)):
        if temperatures[i] >= 273.15:
            inverseTemp.append(1/temperatures[i])
            lnPressure.append(math.log(pressures[i]))
    
    try:
        slope = numpy.polyfit(inverseTemp, lnPressure, 1)[0]
    
        betaGas = -10/slope

    except:
        betaGas = 0
        
    return betaGas

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

def guessPressure(compounds, moleFractions, T):
    noCompounds = len(compounds)
    guessConsts = numpy.array(guessConstants.loc[guessConstants['Compound ID'] == compounds[0]])
    for i in range(noCompounds-1):
        guessConsts = numpy.append(guessConsts, guessConstants.loc[guessConstants['Compound ID'] == compounds[i+1]], axis = 0)

    if T < 273.15:
        guessPressure = 0
        for i in range(noCompounds):
            guessPressure += moleFractions[i]*guessConsts[i][2]*math.exp(guessConsts[i][3]*T)
        
        return guessPressure
    else:
        guessPressure = 0
        for i in range(noCompounds):
            guessPressure += moleFractions[i]*guessConsts[i][4]*math.exp(guessConsts[i][5]*T)

        return guessPressure

def guessTemp(compounds, moleFractions, P):
    noCompounds = len(compounds)
    guessConsts = numpy.array(guessConstants.loc[guessConstants['Compound ID'] == compounds[0]])
    for i in range(noCompounds-1):
        guessConsts = numpy.append(guessConsts, guessConstants.loc[guessConstants['Compound ID'] == compounds[i+1]], axis = 0)

    constantA = 0
    constantB = 0
    for i in range(noCompounds):
        constantA += moleFractions[i]*guessConsts[i][2]
        constantB += moleFractions[i]*guessConsts[i][3]
    
    guessTemp = math.log(P/constantA)/constantB

    if guessTemp >= 273.15:
        constantA = 0
        constantB = 0
        for i in range(noCompounds):
            constantA += moleFractions[i]*guessConsts[i][4]
            constantB += moleFractions[i]*guessConsts[i][5]
    
    guessTemp = math.log(P/constantA)/constantB

    return guessTemp

def hydrationNumber(structure, occupancies):
    if sum(occupancies[0] + occupancies[1]) != 0:
        if structure == "I":
            hydrationNumber = 46/(sum(occupancies[0])*2+sum(occupancies[1])*6)
        else:
            hydrationNumber = 136/(sum(occupancies[0])*16+sum(occupancies[1])*8)
        
        return round(hydrationNumber, 2)
    else:
        return None

def hydrateDensity(structure, occupancies, compoundData, moleFractions, T, P):
    noCompounds = len(moleFractions)
    N_A = 6.022E23
    guestMass = 0
    
    molarmass = numpy.zeros(noCompounds)
    for i in range(noCompounds):
        molarmass[i] = compoundData[i][10]/1000
    
    if structure == "I":
        Vm_water = (11.835+2.217E-5*T+2.242E-6*T**2)**3*(1E-30*N_A/46)-8.006E-9*P/1E6+5.448E-12*(P/1E6)**2 #m3/mol
        waterMass = 18.02/1000/Vm_water #kg/m3
        for i in range(noCompounds):
            guestMass += (molarmass[i]*occupancies[0][i]*2 + molarmass[i]*occupancies[1][i]*6)/Vm_water/46
    elif structure == "II":
        Vm_water = (17.13+2.249E-4*T+2.013E-6*T**2-1.009E-9*T**3)**3*(1E-30*N_A/136)-8.006E-9*P/1E6+5.448E-12*(P/1E6)**2 #m3/mol
        waterMass = 18.02/1000/Vm_water #kg/m3
        for i in range(noCompounds):
            guestMass += (molarmass[i]*occupancies[0][i]*16 + molarmass[i]*occupancies[1][i]*8)/Vm_water/136
            
    return round(waterMass + guestMass, 1)

def equilibriumPressure(temperature, pressure, compounds, moleFractions, saltConcs, inhibitorConcs):
    noCompounds = len(compounds)
    compoundData = numpy.array(fluidProperties.loc[fluidProperties['Compound ID'] == compounds[0]])
    for i in range(noCompounds -1):
        compoundData = numpy.append(compoundData, fluidProperties.loc[fluidProperties['Compound ID'] == compounds[i+1]], axis = 0)

    PvapConsts = numpy.array(vaporPressureConstants.loc[vaporPressureConstants['Compound ID'] == compounds[0]])
    for i in range(noCompounds -1):
        PvapConsts = numpy.append(PvapConsts, vaporPressureConstants.loc[vaporPressureConstants['Compound ID'] == compounds[i+1]], axis = 0)

    kiharaParameters = None

    #Computational Algorithm
    pGuess = pressure
    if temperature > 260 and temperature < 280:
        freezingPoint = 273.15+freezingPointDepression(compounds, temperature, PengRobinson(compoundData, moleFractions, 273.15, pressure)[2], compoundData, pGuess, compoundData[:,7], saltConcs, inhibitorConcs)
        if temperature > freezingPoint:
            waterPhase = "liquid"
        else:
            waterPhase = "ice"
    elif temperature <= 260:
        waterPhase = "ice"
    elif temperature >= 280:
        waterPhase = "liquid"
    
    def f(pressure, temperature):
        vaporFugacities = PengRobinson(compoundData, moleFractions, temperature, abs(pressure))[2]
        f_w = waterFugacity(temperature, pressure, waterPhase, vaporFugacities, compounds, compoundData)
        f_h = hydrateFugacity(temperature, pressure, localPvapConsts, structure, vaporFugacities, compounds, kiharaParameters, compoundData)[0]
        return abs(f_h/f_w-1)
      
    structure = "I"
    pGuess = pressure
    mask = [0 for i in range(len(PvapConsts[:,0]))]
    for i in range(len(PvapConsts[:,0])):
        if sum(numpy.isin(element = PvapConsts[:,1],test_elements=PvapConsts[i,1])) > 1:
            structureMask = numpy.isin(element = PvapConsts[:,0],test_elements="I")
            componentMask = numpy.isin(element = PvapConsts[:,1],test_elements=PvapConsts[i,1])
            mask[i] = numpy.logical_and(structureMask, componentMask)[i]
        else:
            mask[i] = True
        
    localPvapConsts = PvapConsts[mask]

    try:
        if "I" in PvapConsts[:,0]:
            SIEqPressure = abs(scipy.optimize.fsolve(f,pGuess,xtol=errorMargin,args=temperature)[0])
            SIEqFrac = hydrateFugacity(temperature, SIEqPressure, localPvapConsts, structure, PengRobinson(compoundData, moleFractions, temperature, pressure)[2], compounds, kiharaParameters, compoundData)[1]
        else:
            raise
    except:
        SIEqPressure = math.inf
        SIEqFrac = numpy.zeros((2,len(moleFractions)))
    
    pGuess = pressure
    
    structure = "II"
    mask = [0 for i in range(len(PvapConsts[:,0]))]
    for i in range(len(PvapConsts[:,0])):
        if sum(numpy.isin(element = PvapConsts[:,1],test_elements=PvapConsts[i,1])) > 1:
            structureMask = numpy.isin(element = PvapConsts[:,0],test_elements="II")
            componentMask = numpy.isin(element = PvapConsts[:,1],test_elements=PvapConsts[i,1])
            mask[i] = numpy.logical_and(structureMask, componentMask)[i]
        else:
            mask[i] = True
        
    localPvapConsts = PvapConsts[mask]
        
    try:
        SIIEqPressure = abs(scipy.optimize.fsolve(f,[pGuess],xtol=errorMargin,args=temperature)[0])
        SIIEqFrac = hydrateFugacity(temperature, SIIEqPressure, localPvapConsts, structure, PengRobinson(compoundData, moleFractions, temperature, pressure)[2], compounds, kiharaParameters, compoundData)[1]
    except:
        SIIEqPressure = math.inf
        SIIEqFrac = numpy.zeros((2,noCompounds))
    
    if SIIEqPressure >= SIEqPressure:
        eqStructure = "I"
        EqFrac = SIEqFrac
    else:
        eqStructure = "II"
        EqFrac = SIIEqFrac

    eqPressure = min(SIEqPressure, SIIEqPressure)

    if waterPhase == "ice":
        equilPhase = "I-H-V"
    else:
        equilPhase = "L-H-V"

    return eqPressure, eqStructure, EqFrac, hydrationNumber(eqStructure, EqFrac), hydrateDensity(eqStructure, EqFrac, compoundData, moleFractions, temperature, eqPressure), equilPhase

def equilibriumTemperature(temperature, pressure, compounds, moleFractions, saltConcs, inhibitorConcs):
    noCompounds = len(compounds)
    compoundData = numpy.array(fluidProperties.loc[fluidProperties['Compound ID'] == compounds[0]])
    for i in range(noCompounds -1):
        compoundData = numpy.append(compoundData, fluidProperties.loc[fluidProperties['Compound ID'] == compounds[i+1]], axis = 0)

    PvapConsts = numpy.array(vaporPressureConstants.loc[vaporPressureConstants['Compound ID'] == compounds[0]])
    for i in range(noCompounds -1):
        PvapConsts = numpy.append(PvapConsts, vaporPressureConstants.loc[vaporPressureConstants['Compound ID'] == compounds[i+1]], axis = 0)

    kiharaParameters = None

    #Computational Algorithm
    tGuess = temperature
    if temperature > 260 and temperature < 280:
        freezingPoint = 273.15+freezingPointDepression(compounds, tGuess, PengRobinson(compoundData, moleFractions, 273.15, pressure)[2], compoundData, pressure, compoundData[:,7], saltConcs, inhibitorConcs)
        if temperature > freezingPoint:
            waterPhase = "liquid"
        else:
            waterPhase = "ice"
    elif temperature <= 260:
        waterPhase = "ice"
    elif temperature >= 280:
        waterPhase = "liquid"
    
    def f(tGuess, pressure):
        vaporFugacities = PengRobinson(compoundData, moleFractions, abs(tGuess), pressure)[2]
        f_w = waterFugacity(abs(tGuess), pressure, waterPhase, vaporFugacities, compounds, compoundData)
        f_h = hydrateFugacity(abs(tGuess), pressure, localPvapConsts, structure, vaporFugacities, compounds, kiharaParameters, compoundData)[0]
        return abs(f_h/f_w-1)
      
    structure = "I"
    tGuess = temperature
    mask = [0 for i in range(len(PvapConsts[:,0]))]
    for i in range(len(PvapConsts[:,0])):
        if sum(numpy.isin(element = PvapConsts[:,1],test_elements=PvapConsts[i,1])) > 1:
            structureMask = numpy.isin(element = PvapConsts[:,0],test_elements="I")
            componentMask = numpy.isin(element = PvapConsts[:,1],test_elements=PvapConsts[i,1])
            mask[i] = numpy.logical_and(structureMask, componentMask)[i]
        else:
            mask[i] = True
        
    localPvapConsts = PvapConsts[mask]
    
    try:
        if "I" in PvapConsts[:,0]:
            SIEqTemperature = abs(scipy.optimize.fsolve(f,tGuess,xtol=errorMargin,args=pressure)[0])
            SIEqFrac = hydrateFugacity(SIEqTemperature, pressure, localPvapConsts, structure, PengRobinson(compoundData, moleFractions, temperature, pressure)[2], compounds, kiharaParameters, compoundData)[1]
        else:
            raise
    except:
        SIEqTemperature = math.inf
        SIEqFrac = numpy.zeros((2,noCompounds))
    
    tGuess = temperature
    
    structure = "II"
    mask = [0 for i in range(len(PvapConsts[:,0]))]
    for i in range(len(PvapConsts[:,0])):
        if sum(numpy.isin(element = PvapConsts[:,1],test_elements=PvapConsts[i,1])) > 1:
            structureMask = numpy.isin(element = PvapConsts[:,0],test_elements="II")
            componentMask = numpy.isin(element = PvapConsts[:,1],test_elements=PvapConsts[i,1])
            mask[i] = numpy.logical_and(structureMask, componentMask)[i]
        else:
            mask[i] = True
        
    localPvapConsts = PvapConsts[mask]
        
    try:
        SIIEqTemperature = abs(scipy.optimize.fsolve(f,[tGuess],xtol=errorMargin,args=pressure)[0])
        SIIEqFrac = hydrateFugacity(SIIEqTemperature, pressure, localPvapConsts, structure, PengRobinson(compoundData, moleFractions, temperature, pressure)[2], compounds, kiharaParameters, compoundData)[1]
    except:
        SIIEqTemperature = math.inf
        SIIEqFrac = numpy.zeros((2,len(moleFractions)))
    
    if SIIEqTemperature <= SIEqTemperature:
        eqStructure = "I"
        EqFrac = SIEqFrac
    else:
        eqStructure = "II"
        EqFrac = SIIEqFrac

    if SIEqTemperature != math.inf and SIIEqTemperature != math.inf:
        eqTemperature = min(SIEqTemperature, SIIEqTemperature)
    else:
        eqTemperature = math.inf

    if waterPhase == "ice":
        equilPhase = "I-H-V"
    else:
        equilPhase = "L-H-V"

    return eqTemperature, eqStructure, EqFrac, hydrationNumber(eqStructure, EqFrac), hydrateDensity(eqStructure, EqFrac, compoundData, moleFractions, eqTemperature, pressure), equilPhase