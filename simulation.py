import math 
import matplotlib.pyplot as plt
import numpy
import pandas
import scipy
import warnings

warnings.filterwarnings('ignore')

#Extract data from excel files
fluidProperties = pandas.read_excel('Data.xlsx', sheet_name='Fluid Properties')
hydrateCellProperties = pandas.read_excel('Data.xlsx', sheet_name='Hydrate Cell Properties')
kiharaCellParameters = pandas.read_excel('Data.xlsx', sheet_name='Kihara Cell Parameters')
vaporPressureConstants = pandas.read_excel('Data.xlsx', sheet_name='Vapor Pressure Constants')
mixConstants = pandas.read_excel('Data.xlsx', sheet_name='Binary Interaction Parameters')

#Constants
errorMargin = 1E-4
R = 8.31446261815324 #m^3 Pa/mol K
boltzmannConstant = 1.380649E-23 #m^2 kg/s^2 K

#PRSV Equation of State (DONE)
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
        kappa = 0.37464 + 1.54226*omega - 0.26992*(omega**2)
        alpha = (1 + kappa*(1-math.sqrt(T/Tc)))**2
        a[i] = 0.45724*((R**2)*(Tc**2)/Pc)*alpha
        
    #Obtain the interaction parameters for all combinations of compounds
    interactionParameters = [[0 for i in range(len(moleFractions))] for j in range(len(moleFractions))]
    for i in range(len(moleFractions)):
        for j in range(len(moleFractions)):
            interactionParameters[i][j] = mixConstants.loc[(mixConstants['Compound 1'] == compoundData[i][1]), (compoundData[j][1])].reset_index(drop=True)[0]
        
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
        
    return VmVap, VmLiq, fugVap#, fugLiq

#For a given hydrate structure, return its radii and coordination numbers (DONE)
def getHydrateCellProperties(structure):
    cellProperties = numpy.array(hydrateCellProperties.loc[(hydrateCellProperties['Structure'] == structure)])
    return cellProperties
   
#Equation 7 (DONE)
def delta(N, r, RCell, a):
    delta = ((1-r/RCell-a/RCell)**(-1*N)-(1+r/RCell-a/RCell)**(-1*N))/N
    return delta 
   
#Equation 6 (DONE)
def W(r, RCell, z, epsilon, sigma, a):
    d10 = delta(10, r, RCell, a)
    d11 = delta(11, r, RCell, a)
    d4 = delta(4, r, RCell, a)
    d5 = delta(5, r, RCell, a)
    W = 2*z*epsilon*((sigma**12)/(r*RCell**11)*(d10+a/RCell*d11)-sigma**6/(r*RCell**5)*(d4+a/RCell*d5))
    return W  
   
#Equation 4 (DONE)
def Lang_Const(T, cellRadii, a, RCell, z, epsilon, sigma):
    def integrand(r):
        W1 = W(r, cellRadii[0], z[0], epsilon, sigma, a)
        W2 = W(r, cellRadii[1], z[1], epsilon, sigma, a)
        W3 = W(r, cellRadii[2], z[2], epsilon, sigma, a)
        x = math.exp(-1*(W1+W2+W3)/T)*r**2
        return x
    
    Cml = 4*math.pi/(boltzmannConstant*T)*(scipy.integrate.quad(integrand, 0, RCell-2*a)[0])
    return Cml

#Equation 3 #TODO: FIX THIS
def frac(T, kiharaParameters, cellProperties, vaporFugacityCoeffs):
    fractions = [[0 for i in range(len(vaporFugacityCoeffs))], [0 for i in range(len(vaporFugacityCoeffs))]]
    denominator = 0 #should this be reset for each cage?
    for i in range(2):
        
        #Hydrate Cell Properties
        cellRadii = [0, 0, 0]
        for j in range(3):
            cellRadii[j] = cellProperties[i][2+j]*1E-10
        z = [0, 0, 0]
        for j in range(3):
            z[j] = cellProperties[i][5+j]
        RCell = cellProperties[i][8]*1E-10
               
        langConsts = [0 for j in range(len(vaporFugacityCoeffs))]
        for j in range(len(vaporFugacityCoeffs)):
            #Kihara Cell Parameters
            epsilon = math.sqrt(102.134*kiharaParameters[j][2])
            sigma = (kiharaParameters[j][3]+ 3.56438)/2*1E-10
            a = kiharaParameters[j][4]/2*1E-10
            langConsts[j] = Lang_Const(T, cellRadii, a, RCell, z, epsilon, sigma)
            denominator += langConsts[j]*vaporFugacityCoeffs[j]
            
        for j in range(len(vaporFugacityCoeffs)):    
            fractions[i][j] = langConsts[j]*vaporFugacityCoeffs[j]/(1 + denominator)
    
    return fractions

#Equation 2 (DONE)
def deltaHydratePotential(T, kiharaParameters, structure, vaporFugacityCoeffs):
    cellProperties = getHydrateCellProperties(structure)
    fractions = 0
    Deltamu_H_w = 0
    
    fractions = frac(T, kiharaParameters, cellProperties, vaporFugacityCoeffs)
    Deltamu_H_w += -1*R*T*(cellProperties[0][10])*math.log(1-sum(fractions[0]))
    Deltamu_H_w += -1*R*T*(cellProperties[1][10])*math.log(1-sum(fractions[1]))
    return Deltamu_H_w

#Equation 18 (DONE)
def henrysLawConst(compound, T):
    henrysLawConstants = pandas.read_excel('Data.xlsx', sheet_name='Henrys Law Parameters')
    constants = numpy.array(henrysLawConstants.loc[henrysLawConstants['Compound ID'] == compound])[0]
    H_i = 101325*math.exp(-1*(constants[2] + constants[3]/T + constants[4]*math.log(T) + constants[5]*T))
    return H_i

#Equation 17 (DONE)
def liqPhaseComposition(compounds, T, fug_vap, Vi, P, Psat):
    x = [0 for i in range(len(compounds)+1)]
    for i in range(len(fug_vap)):
        x[i+1] = fug_vap[i]/(henrysLawConst(compounds[i], T)*math.exp(Vi[i]*(P-Psat)/(R*T)))
    x[0] = 1-sum(x) #Water composition
    return x

#Equation 19
def freezingPointDepression():
    deltadT = R*(273.15)**2/6011*math.log(1)
    return deltadT

def activityCoeff():
    return 1

#Equations 15 and 16
def waterFugacity(T, P, phase, fug_vap, compounds, Vi):
    if phase == "ice":
        Vm_water = 1.912E-5 + T*8.387E-10 + (T**2)*4.016E-12
        Psat_water = math.exp(4.6056*math.log(T)-5501.1243/T+2.9446-T*8.1431E-3)
        f_w = Psat_water*math.exp(Vm_water*(P-Psat_water)/(R*T))
    elif phase == "liquid":
        Vm_water = math.exp(-10.921 + 2.5E-4*(T-273.15) - 3.532E-4*(P/1E6-0.101325) + 1.559E-7*((P/1E6-.101325)**2))
        Psat_water = math.exp(4.1539*math.log(T)-5500.9332/T+7.6537-16.1277E-3*T)
        f_w = liqPhaseComposition(compounds, T, fug_vap, Vi, P, Psat_water)[0]*activityCoeff()*Psat_water*math.exp(Vm_water*(P-Psat_water)/(R*T))
    return f_w

#Equation A2
def hydrateFugacity(T, P, PvapConsts, structure, fug_vap, compounds, kiharaParameters):
    N_A = 6.022E23
    if structure == "I":
        Vm_water = (11.835+2.217E-5*T+2.242E-6*T**2)*(1E-30*N_A/46)-8.006E-9*P/1E6+5.448E-12*(P/1E6)**2
        #Vm_water = (11.835+2.217E-5*T+2.242E-6*T**2)**3*(1E-30*N_A/46)
    elif structure == "II":
        Vm_water = (17.13+2.249E-4*T+2.013E-6*T**2+1.009E-9*T**3)*(1E-30*N_A/136)-8.006E-9*P/1E6+5.448E-12*(P/1E6)**2
    A = PvapConsts[0][3]
    B = PvapConsts[0][4] #Figure out how this works with mixing
    D = PvapConsts[0][6]
    Psat_water = math.exp(A*math.log(T)+B/T+2.7789+D*T)
    f_h = Psat_water*math.exp(Vm_water*(P-Psat_water)/(R*T))*math.exp(-1*deltaHydratePotential(T, kiharaParameters, structure, fug_vap)/(R*T))
    return f_h

def equilibriumPressure(temperature, pressure, compounds, moleFractions):
    compoundData = numpy.array(fluidProperties.loc[fluidProperties['Compound ID'] == compounds[0]])
    for i in range(len(compounds)-1):
        compoundData = numpy.append(compoundData, fluidProperties.loc[fluidProperties['Compound ID'] == compounds[i+1]], axis = 0)

    PvapConsts = numpy.array(vaporPressureConstants.loc[vaporPressureConstants['Compound ID'] == compounds[0]])
    for i in range(len(compounds)-1):
        PvapConsts = numpy.append(vaporPressureConstants, vaporPressureConstants.loc[vaporPressureConstants['Compound ID'] == compounds[i+1]], axis = 0)

    kiharaParameters = numpy.array(kiharaCellParameters.loc[kiharaCellParameters['Compound ID'] == compounds[0]])
    for i in range(len(compounds)-1):
        kiharaParameters = numpy.append(kiharaParameters, kiharaCellParameters.loc[kiharaCellParameters['Compound ID'] == compounds[i+1]], axis = 0)

    #Computational Algorithm
    structure = "I"
    pGuess = pressure

    if temperature > 273.15-freezingPointDepression():
        waterPhase = "liquid"
    else:
        waterPhase = "ice"
        
    eosResult = PengRobinson(compoundData, moleFractions, temperature, pressure)
    vaporFugacities = eosResult[2]
        
    f_w = waterFugacity(temperature, pressure, waterPhase, vaporFugacities, compounds, compoundData[:,6])
        
    f_h = hydrateFugacity(temperature, pressure, PvapConsts, structure, vaporFugacities, compounds, kiharaParameters)

    while abs(f_h/f_w - 1) > errorMargin:
        pGuess = pGuess*(f_h/f_w)
        
        vaporFugacities = PengRobinson(compoundData, moleFractions, temperature, pGuess)[2]
        f_w = waterFugacity(temperature, pGuess, waterPhase, vaporFugacities, compounds, compoundData[:,6])
        f_h = hydrateFugacity(temperature, pGuess, PvapConsts, structure, vaporFugacities, compounds, kiharaParameters)

    SIEqPressure = pGuess

    structure = "II"
    pGuess = pressure

    eosResult = PengRobinson(compoundData, moleFractions, temperature, pressure)
    vaporFugacities = eosResult[2]
        
    f_w = waterFugacity(temperature, pressure, waterPhase, vaporFugacities, compounds, compoundData[:,6])
        
    f_h = hydrateFugacity(temperature, pressure, PvapConsts, structure, vaporFugacities, compounds, kiharaParameters)

    while abs(f_h/f_w - 1) > errorMargin:
        pGuess = pGuess*(f_h/f_w)
            
        vaporFugacities = PengRobinson(compoundData, moleFractions, temperature, pGuess)[2]
        f_w = waterFugacity(temperature, pGuess, waterPhase, vaporFugacities, compounds, compoundData[:,6])
        f_h = hydrateFugacity(temperature, pGuess, PvapConsts, structure, vaporFugacities, compounds, kiharaParameters)

    SIIEqPressure = pGuess
    
    return min(SIEqPressure, SIIEqPressure)

#Runtime-----------------------------------------------------------------------
minTemp = float(input("Lower Temperature (K): ")) #Temp in K
maxTemp = float(input("Upper Temperature (K): "))
pressure = float(input("Pressure (MPa): "))*1E6 #Pressure in Pa
numberOfCompounds = int(input("No. of Compounds Excluding Water: "))
compounds = []
moleFractions = []
print("\n 1. Methane      2. Ethane\n 3. Propane      4. i-Butane\n 5. c-C3H6       6. H2S\n 7. Nitrogen     8. CO2\n")
for i in range(numberOfCompounds):
    compounds += [int(input("Compound ID " + str(i + 1) + " : "))]
    if numberOfCompounds > 1:
        moleFractions += [float(input("Mole Fraction of Compound " + str(i + 1) + " : "))]
    else:
        moleFractions = [1]
noPoints = int(input("Number of Data Points: "))
    
T = numpy.arange(minTemp, maxTemp+(maxTemp-minTemp)/noPoints, (maxTemp-minTemp)/(noPoints-1))
eqPressure = [0 for i in range(len(T))]
for i in range(len(T)):
    eqPressure[i] = equilibriumPressure(T[i], pressure, compounds, moleFractions)/1E6 #In MPa

plt.plot(T, eqPressure, '-ok')
plt.yscale("log")
plt.title("Equilibrium Predictions")
plt.xlabel("Temperature (K)")
plt.ylabel("Pressure (MPa)")
plt.show