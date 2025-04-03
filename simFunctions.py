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
hydrateCellProperties = pandas.read_excel('Data.xlsx', sheet_name='Hydrate Cell Properties')
#kiharaCellParameters = pandas.read_excel('Data.xlsx', sheet_name='Kihara Cell Parameters')
vaporPressureConstants = pandas.read_excel('Data.xlsx', sheet_name='Vapor Pressure Constants')
mixConstants = pandas.read_excel('Data.xlsx', sheet_name='Binary Interaction Parameters')
saltData = pandas.read_excel('Data.xlsx', sheet_name='Salts Data').to_numpy()
inhibitorData = pandas.read_excel('Data.xlsx', sheet_name='Inhibitor Data').to_numpy()

#Constants
errorMargin = 1E-9
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
    #Calculate "b" value for each pure substance and add to weighted and unweighted lists
    bmix = 0
    b = [0]*len(moleFractions)
    for i in range(len(moleFractions)):    
        Tc = compoundData[i][2]
        Pc = compoundData[i][3]*1E6
        b[i] = 0.07780*(R*Tc/Pc)
        bmix += b[i]*moleFractions[i]
        
    #Calculate the "a" values for all pure compounds
    a = [0]*len(moleFractions)
    for i in range(len(moleFractions)):
        Tc = compoundData[i][2]
        Pc = compoundData[i][3]*1E6
        omega = compoundData[i][4]
        kappa1 = compoundData[i][5]
        kappa = 0.378893 + 1.4897153*omega - 0.17131848*omega**2 + 0.0196554*omega**3 + kappa1*(1+math.sqrt(T/Tc))*(0.7-(T/Tc))
        #kappa = 0.37464 + 1.54226*omega - 0.26992*(omega**2)
        alpha = (1 + kappa*(1-math.sqrt(T/Tc)))**2
        a[i] = 0.45724*((R**2)*(Tc**2)/Pc)*alpha
        
    #Obtain the interaction parameters for all combinations of compounds
    if len(moleFractions) > 1:
        interactionParameters = [[0 for i in range(len(moleFractions))] for j in range(len(moleFractions))]
        for i in range(len(moleFractions)):
            for j in range(len(moleFractions)):
                interactionParameters[i][j] = mixConstants.loc[(mixConstants['Compound 1'] == compoundData[i][1]), (compoundData[j][1])].reset_index(drop=True)[0]
    else:
        interactionParameters = [[0]]
        
    #Finally, determine the partial "a" values based on mole ratios
    amix = 0
    xia = [0 for i in range(len(moleFractions))]
    for i in range(len(moleFractions)):
        for j in range(len(moleFractions)):
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
    
    fugVap = [0 for i in range(len(moleFractions))]
    #fugLiq = [0 for i in range(len(compounds))]
    
    #Calculate fugacities for each of the compounds
    for i in range(len(moleFractions)):
        fugVap[i] = P*moleFractions[i]*math.exp((b[i]/bmix)*(ZVmix-1)-math.log(ZVmix-B)-(A/(2*math.sqrt(2)*B))*(2*xia[i]/amix-b[i]/bmix)*math.log((ZVmix+(1+math.sqrt(2))*B)/(ZVmix+(1-math.sqrt(2))*B)))
        #fugLiq[i] = P*moleFractions[i]*math.exp((b[i]/bmix)*(ZLmix-1)-math.log(ZLmix-B)-(A/(2*math.sqrt(2)*B))*(2*xia[i]/amix-b[i]/bmix)*math.log((ZLmix+(1+math.sqrt(2))*B)/(ZLmix+(1-math.sqrt(2))*B)))
       
    return VmVap, VmLiq, fugVap, ZVmix

#For a given hydrate structure, return its radii and coordination numbers
def getHydrateCellProperties(structure):
    cellProperties = numpy.array(hydrateCellProperties.loc[(hydrateCellProperties['Structure'] == structure)])
    return cellProperties
   
def delta(N, r, RCell, a):
    delta = ((1-r/RCell-a/RCell)**(-1*N)-(1+r/RCell-a/RCell)**(-1*N))/N
    return delta 
   
def W(r, RCell, z, epsilon, sigma, a):
    d10 = delta(10, r, RCell, a)
    d11 = delta(11, r, RCell, a)
    d4 = delta(4, r, RCell, a)
    d5 = delta(5, r, RCell, a)
    W = 2*z*epsilon*((sigma**12)/(r*RCell**11)*(d10+a/RCell*d11)-sigma**6/(r*RCell**5)*(d4+a/RCell*d5))
    return W  

#Calculates Langmuir Constant Cml
def Lang_Const(T, cellRadii, a, RCell, z, epsilon, sigma, structure, shell, compound):
    '''def integrand(r):
        W1 = W(r, cellRadii[0], z[0], epsilon, sigma, a)
        W2 = W(1E-14, cellRadii[1], z[1], epsilon, sigma, a)
        W3 = W(1E-14, cellRadii[2], z[2], epsilon, sigma, a)
        x = math.exp(-1*(W1+W2+W3)/T)*(r**2)
        return x
    
    Cml = 4*math.pi/(boltzmannConstant*T)*(scipy.integrate.quad(integrand, 0, RCell-a)[0])'''
    
    langParameters = pandas.read_excel('Data.xlsx', sheet_name='Langmuir Parameters')
    
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
    Cgg = [[1,1] for i in range(len(I))]
    
    C6 = [[0 for i in range(len(I))] for i in range(len(I))]
    C8 = [[0 for i in range(len(I))] for i in range(len(I))]
    C10 = [[0 for i in range(len(I))] for i in range(len(I))]
            
    for i in range(len(I)):
        for j in range(len(I)):
            C6[i][j] = (3/2)*alpha[i]*alpha[j]*I[i]*I[j]/(I[i]+I[j])*23.05/.001987
            C8[i][j] = 496778.3824*alpha[i]*alpha[j]*(I[i]/(2*I[i]+I[j])+I[j]/(2*I[j]+I[i]))
            C10[i][j] = 13260598.42*alpha[i]*alpha[j]/(I[i]+I[j])
    
    integrationConstants = pandas.read_excel('Data.xlsx', sheet_name='A_int')
    
    wrgg = [[[[0,0],[0,0]] for i in range(len(I))] for i in range(len(I))]
    
    #Unsure exactly what these are relative to r
    if structure == "I":
        a_lc = 12.03
    elif structure == "II":
        a_lc = 17.31
        
    def Aij(i, j, k, l):
        condition = (integrationConstants['Structure'] == i) & \
            (integrationConstants['Center Cage'] == j) & \
            (integrationConstants['Parameter'] == k) & \
            (integrationConstants['Like/Dislike'] == l)
        Aij = integrationConstants.loc[condition, "A"].values[0]
        return Aij
        
    for i in range(len(I)):
        for j in range(len(I)):
            for k in range(2):
                r0 = a_lc*Aij(structure, k, 3, 0)
                wrgg[i][j][k][0] = -1*C6[i][j]*Aij(structure, k, 0, 0)/(r0**6)-C8[i][j]*Aij(structure, k, 1, 0)/(r0**8)-C10[i][j]*Aij(structure, k, 2, 0)/(r0**10)
                r1 = a_lc*Aij(structure, k, 3, 1)
                wrgg[i][j][k][1] = -1*C6[i][j]*Aij(structure, k, 0, 1)/(r1**6)-C8[i][j]*Aij(structure, k, 1, 1)/(r1**8)-C10[i][j]*Aij(structure, k, 2, 1)/(r1**10)
    
    #Obtain the interaction energy of guest
    for i in range(len(I)):
        for j in range(len(I)):
            for k in range(2):
                Cgg[i][k] *= math.exp(-1*(wrgg[i][j][k][0]*fracs[0][j])/T)*math.exp(-1*(wrgg[i][j][k][1]*fracs[1][j])/T)
    
    return Cgg

#Calculates the fractional occupancy of small and large shells by each component
def frac(T, kiharaParameters, vaporFugacities, compoundData, structure, compounds):
    cellProperties = getHydrateCellProperties(structure) 
    guessFractions = [[0 for i in range(len(vaporFugacities))],[0 for i in range(len(vaporFugacities))]]
    oldGuessFractions = [[0 for i in range(len(vaporFugacities))],[0 for i in range(len(vaporFugacities))]]
    Cgg = [[1, 1] for i in range(len(vaporFugacities))]
    langConsts = [[None, None] for i in range(len(vaporFugacities))]
    RCell = [0,0]
    a = [0 for i in range(len(vaporFugacities))]
    fracDiff = 0.5
    while abs(fracDiff) >= 1E-4:
        fracDiff = 0
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
                   
            for j in range(len(vaporFugacities)):
                #Kihara Cell Parameters
                epsilon = None #math.sqrt(102.134*kiharaParameters[j][2])
                sigma = None #(kiharaParameters[j][3]+ 3.56438)/2*1E-10
                a[j] = None #kiharaParameters[j][4]/2*1E-10
                langConsts[j][i] = Lang_Const(T, cellRadii, a[j], RCell[i], z, epsilon, sigma, structure, i, compounds[j])*(1-fracDiff)
                denominator += Cgg[j][i]*langConsts[j][i]*vaporFugacities[j]
                    
            for j in range(len(vaporFugacities)):
                guessFractions[i][j] = Cgg[j][i]*langConsts[j][i]*vaporFugacities[j]/(1 + denominator)
                fracDiff += oldGuessFractions[i][j]-guessFractions[i][j]
                oldGuessFractions[i][j]=guessFractions[i][j]
            
        Cgg = Lang_GG_Const(T, compoundData, guessFractions, structure)

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
    henrysLawConstants = pandas.read_excel('Data.xlsx', sheet_name='Henrys Law Parameters')
    constants = numpy.array(henrysLawConstants.loc[henrysLawConstants['Compound ID'] == compound])[0]
    H_i = 101325*math.exp(-1*(constants[2]/1.987 + constants[3]/T/1.987 + constants[4]*math.log(T)/1.987 + constants[5]*T/1.987))
    return H_i

#Infinite Dilution Compressibility Factor
def Z(compoundData, T, P):
    Z = [0 for i in range(len(compoundData))]
    waterData = numpy.array(fluidProperties.loc[fluidProperties['Compound ID'] == 0])[0]
    for i in range(len(compoundData)):
        localCompoundData = compoundData[i]
        localCompoundData = numpy.column_stack((localCompoundData, waterData))
        localCompoundData = localCompoundData.transpose()
        Z[i] = PengRobinson(localCompoundData, [0.00001, 0.99999], 273.15, P)[3]
    return Z

#Equation 17
def liqPhaseComposition(compounds, T, fug_vap, compoundData, P, Psat):
    x = [0 for i in range(len(compounds)+1)]
    for i in range(len(fug_vap)):
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
    for i in range(len(compounds)):
        denominator += frac[0][i]*cellProperties[0][10]
        denominator += frac[1][i]*cellProperties[1][10]
            
    Z = [0 for i in range(len(compounds))]
    for i in range(len(compounds)):
        for j in range(2):
            Z[i] += (frac[j][i]*cellProperties[j][10])/denominator
            
    for i in range(len(compounds)):
        A += PvapConsts[i,3]*Z[i]
        B += PvapConsts[i,4]*Z[i]
        D += PvapConsts[i,6]*Z[i]
    
    Psat_water = math.exp(A*math.log(T)+B/T+2.7789+D*T)
    
    f_h = Psat_water*math.exp(Vm_water*(P-Psat_water)/(R*T))*math.exp(dH[0])
    return f_h,frac

def HuLeeSum(T, saltConcs, inhibitorConcs, betaGas):
    if T >= 273.15:
        if sum(saltConcs) != 0:
            catMols = [0 for i in range(len(saltConcs))]
            anMols = [0 for i in range(len(saltConcs))]
            catMolFracs = [0 for i in range(len(saltConcs))]
            anMolFracs = [0 for i in range(len(saltConcs))]
            for i in range(len(saltConcs)):
                catMols[i] = saltConcs[i]/saltData[i][1]*saltData[i][4]
                anMols[i] = saltConcs[i]/saltData[i][1]*saltData[i][5]
            totalSaltMols = sum(catMols)+sum(anMols)
            waterMols = (100-sum(saltConcs))/18.015
            for i in range(len(saltConcs)):
                catMolFracs[i] = saltData[i][2]*catMols[i]/(waterMols+totalSaltMols)
                anMolFracs[i] = saltData[i][3]*anMols[i]/(waterMols+totalSaltMols)
                
            inhibitorMols = [0 for i in range(len(inhibitorConcs))]
            inhibitorMolFracs = [0 for i in range(len(inhibitorConcs))]
            for i in range(len(inhibitorConcs)):
                inhibitorMols[i] = inhibitorConcs[i]/(1-0.01*inhibitorConcs[i])/inhibitorData[i][1]
            totalMols = totalSaltMols + sum(inhibitorMols)
            for i in range(len(inhibitorConcs)):
                inhibitorMolFracs[i] = inhibitorMols[i]/totalMols

            totalEffSaltMols = sum(catMolFracs) + sum(anMolFracs)
        else:
            waterMols = (100-sum(saltConcs))/18.015
            totalSaltMols = 0
            totalEffSaltMols = 0

        inhibitorMols = [0 for i in range(len(inhibitorConcs))]
        inhibitorMolFracs = [0 for i in range(len(inhibitorConcs))]
        for i in range(len(inhibitorConcs)):
            inhibitorMols[i] = inhibitorConcs[i]/(1-0.01*inhibitorConcs[i])/inhibitorData[i][1]
        totalMols = totalSaltMols + sum(inhibitorMols)
        for i in range(len(inhibitorConcs)):
            inhibitorMolFracs[i] = inhibitorMols[i]/(totalMols+waterMols)

        lnawSalts = -1.06152*totalEffSaltMols+3.25726*totalEffSaltMols*totalEffSaltMols-37.2263*totalEffSaltMols*totalEffSaltMols*totalEffSaltMols
        lnawInhibitors = inhibitorMols = [0 for i in range(len(inhibitorConcs))]
        for i in range(len(inhibitorConcs)):
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
    inhibitorConcs = [0 for i in range(noInhibitors)]
    saltConcs = [0 for i in range(noSalts)]
    
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

def tempConversion(tempUnit, T, isInverted):
    if isInverted == False:
        if tempUnit == "K":
            T = T
        elif tempUnit == "°C":
            T = T - 273.15
        elif tempUnit == "°F":
            T = (T-32)/1.8 + 273.15
        elif tempUnit == "R":
            T = T/1.8
    else:
        if tempUnit == "K":
            T = T
        elif tempUnit == "°C":
            T = T + 273.15
        elif tempUnit == "°F":
            T = (T-273.15)*1.8 + 32
        elif tempUnit == "R":
            T = T*1.8
    return T
    
def pressureConversion(pressureUnit, P, isInverted):
    if isInverted == False:
        if pressureUnit == "MPa":
            P = P
        elif pressureUnit == "psia":
            P /= 145.038
        elif pressureUnit == "bar":
            P /= 10
        elif pressureUnit == "mmHg":
            P /= 7500.62
        #P *= 1E6 #MPa to Pa
    else:
        P /= 1E6 #Pa to MPa   
        if pressureUnit == "MPa":
            P = P
        elif pressureUnit == "psia":
            P *= 145.038
        elif pressureUnit == "bar":
            P *= 10
        elif pressureUnit == "mmHg":
            P *= 7500.62
    return P

def guessPressure(compounds, moleFractions, T):
    compoundData = numpy.array(fluidProperties.loc[fluidProperties['Compound ID'] == compounds[0]])
    for i in range(len(compounds)-1):
        compoundData = numpy.append(compoundData, fluidProperties.loc[fluidProperties['Compound ID'] == compounds[i+1]], axis = 0)

    pressureConstant = 0
    for i in range(len(compoundData)):
        pressureConstant += moleFractions[i]*compoundData[i][6]
    
    guessPressure = math.exp(pressureConstant*T)
    return guessPressure

def guessTemp(compounds, moleFractions, P):
    compoundData = numpy.array(fluidProperties.loc[fluidProperties['Compound ID'] == compounds[0]])
    for i in range(len(compounds)-1):
        compoundData = numpy.append(compoundData, fluidProperties.loc[fluidProperties['Compound ID'] == compounds[i+1]], axis = 0)

    pressureConstant = 0
    for i in range(len(compoundData)):
        pressureConstant += moleFractions[i]*compoundData[i][6]
    
    guessTemp = math.log(P)/pressureConstant
    return guessTemp

def hydrationNumber(structure, occupancies):
    if structure == "I":
        hydrationNumber = 46/(sum(occupancies[0])*2+sum(occupancies[1])*6)
    else:
        hydrationNumber = 136/(sum(occupancies[0])*16+sum(occupancies[1])*8)
    
    return round(hydrationNumber, 2)

def hydrateDensity(structure, occupancies, compoundData, moleFractions, T, P):
    N_A = 6.022E23
    guestMass = 0
    
    molarmass = [0 for i in range(len(moleFractions))]
    for i in range(len(moleFractions)):
        molarmass[i] = compoundData[i][10]/1000
    
    if structure == "I":
        Vm_water = (11.835+2.217E-5*T+2.242E-6*T**2)**3*(1E-30*N_A/46)-8.006E-9*P/1E6+5.448E-12*(P/1E6)**2 #m3/mol
        waterMass = 18.02/1000/Vm_water #kg/m3
        for i in range(len(moleFractions)):
            guestMass += (molarmass[i]*occupancies[0][i]*2 + molarmass[i]*occupancies[1][i]*6)/Vm_water/46
    elif structure == "II":
        Vm_water = (17.13+2.249E-4*T+2.013E-6*T**2-1.009E-9*T**3)**3*(1E-30*N_A/136)-8.006E-9*P/1E6+5.448E-12*(P/1E6)**2 #m3/mol
        waterMass = 18.02/1000/Vm_water #kg/m3
        for i in range(len(moleFractions)):
            guestMass += (molarmass[i]*occupancies[0][i]*16 + molarmass[i]*occupancies[1][i]*8)/Vm_water/136
            
    return round(waterMass + guestMass, 2)

def equilibriumPressure(temperature, pressure, compounds, moleFractions, saltConcs, inhibitorConcs):
    compoundData = numpy.array(fluidProperties.loc[fluidProperties['Compound ID'] == compounds[0]])
    for i in range(len(compounds)-1):
        compoundData = numpy.append(compoundData, fluidProperties.loc[fluidProperties['Compound ID'] == compounds[i+1]], axis = 0)

    PvapConsts = numpy.array(vaporPressureConstants.loc[vaporPressureConstants['Compound ID'] == compounds[0]])
    for i in range(len(compounds)-1):
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
        return abs(f_h-f_w)
      
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
        SIEqFrac = [[0 for i in range(len(moleFractions))],[0 for i in range(len(moleFractions))]]
    
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
        if "II" in PvapConsts[:,0]:
            SIIEqPressure = abs(scipy.optimize.fsolve(f,[pGuess],xtol=errorMargin,args=temperature)[0])
            SIIEqFrac = hydrateFugacity(temperature, SIIEqPressure, localPvapConsts, structure, PengRobinson(compoundData, moleFractions, temperature, pressure)[2], compounds, kiharaParameters, compoundData)[1]
        else:
            raise
    except:
        SIIEqPressure = math.inf
        SIIEqFrac = [[0 for i in range(len(moleFractions))],[0 for i in range(len(moleFractions))]]
    
    if SIIEqPressure >= SIEqPressure:
        eqStructure = "I"
        EqFrac = SIEqFrac
    else:
        eqStructure = "II"
        EqFrac = SIIEqFrac

    eqPressure = min(SIEqPressure, SIIEqPressure)

    return eqPressure, eqStructure, EqFrac, hydrationNumber(eqStructure, EqFrac), hydrateDensity(eqStructure, EqFrac, compoundData, moleFractions, temperature, eqPressure)

def equilibriumTemperature(temperature, pressure, compounds, moleFractions, saltConcs, inhibitorConcs):
    compoundData = numpy.array(fluidProperties.loc[fluidProperties['Compound ID'] == compounds[0]])
    for i in range(len(compounds)-1):
        compoundData = numpy.append(compoundData, fluidProperties.loc[fluidProperties['Compound ID'] == compounds[i+1]], axis = 0)

    PvapConsts = numpy.array(vaporPressureConstants.loc[vaporPressureConstants['Compound ID'] == compounds[0]])
    for i in range(len(compounds)-1):
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
        f_w = waterFugacity(tGuess, pressure, waterPhase, vaporFugacities, compounds, compoundData)
        f_h = hydrateFugacity(tGuess, pressure, localPvapConsts, structure, vaporFugacities, compounds, kiharaParameters, compoundData)[0]
        return abs(f_h-f_w)
      
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
        SIEqFrac = [[0 for i in range(len(moleFractions))],[0 for i in range(len(moleFractions))]]
    
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
        if "II" in PvapConsts[:,0]:
            SIIEqTemperature = abs(scipy.optimize.fsolve(f,[tGuess],xtol=errorMargin,args=pressure)[0])
            SIIEqFrac = hydrateFugacity(SIIEqTemperature, pressure, localPvapConsts, structure, PengRobinson(compoundData, moleFractions, temperature, pressure)[2], compounds, kiharaParameters, compoundData)[1]
        else:
            raise
    except:
        SIIEqTemperature = math.inf
        SIIEqFrac = [[0 for i in range(len(moleFractions))],[0 for i in range(len(moleFractions))]]
    
    if SIIEqTemperature <= SIEqTemperature:
        eqStructure = "I"
        EqFrac = SIEqFrac
    else:
        eqStructure = "II"
        EqFrac = SIIEqFrac

    eqTemperature = min(SIEqTemperature, SIIEqTemperature)

    return eqTemperature, eqStructure, EqFrac, hydrationNumber(eqStructure, EqFrac), hydrateDensity(eqStructure, EqFrac, compoundData, moleFractions, eqTemperature, pressure)