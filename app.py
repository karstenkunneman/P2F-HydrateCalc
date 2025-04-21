#This work is licensed under the Creative Commons Attribution-ShareAlike 4.0 International License. To view a copy of this license, visit https://creativecommons.org/licenses/by-sa/4.0/.
#Created by Karsten Kunneman and Amadeu K. Sum at the Colorado School of Mines
#©2025, All Rights Reserved

import simFunctions
import streamlit as st
import pandas as pd
import numpy
import matplotlib.pyplot as plt
import time
import math

IDs, compounds = simFunctions.getComponents()
componentList = []
for i in range(len(IDs)):
    componentList.append(compounds[i])

st.title('Phases to Flow :: Gas Hydrate Equilibrium Prediction Calculator')
st.caption('Version 2025-04-16')

programType = st.radio("Calculation Type", ["Equilibrium Calculation", "Minimum Concentration Calculation"], horizontal=True)

if programType == "Equilibrium Calculation":
    #Mole Fraction Input Table
    components = []
    moleFractions = []

    compDf = pd.DataFrame([])
    compDf = pd.concat([compDf, pd.DataFrame([{'Component': componentList[0], 'Mole Fraction': 1.}])], ignore_index=True)
    for i in range(len(componentList)-1):
        compDf = pd.concat([compDf, pd.DataFrame([{'Component': componentList[i+1], 'Mole Fraction': 0.}])], ignore_index=True)

    inputCompDf = st.data_editor(compDf, hide_index=True)

    moleFracInput = inputCompDf['Mole Fraction'].tolist()

    for i in range(len(componentList)):
        if moleFracInput[i] > 0:
            components.append(i + 1)
            moleFractions.append(moleFracInput[i])

    #User temperature and guess pressure input file w/ template
    csvGuesses = st.file_uploader("Upload Temperatures and Guess Pressures (Optional)", ['csv'])
    csvTemplate = st.download_button("Guess File Template", open("Input Template.csv", encoding='utf-8'), file_name="Input Template.csv")

    #Unit Selector
    tempUnit = st.radio("Temperature Unit", ["K", "°C", "°F"], horizontal=True)
    pressureUnit = st.radio("Pressure Unit", ["MPa", "bar", "psia"], horizontal=True)

    #Defined Variable Selector
    definedVariable = st.radio("Defined Variable", ["T", "P"], horizontal=True)

    #If no input file given, take inputs from user through UI
    if csvGuesses == None:
        calculateRange = st.toggle('Calculation Range', False)
        userGuess = st.toggle('Manually Guess?', False)
        if calculateRange == True:
            if definedVariable == "T":
                c1, c2 = st.columns(2)
                with c1:
                    minTemp = float(st.text_input('Minimum Temperature ('+tempUnit+')', round(simFunctions.tempConversion(tempUnit, 240.2, False), 1)))
                with c2:
                    maxTemp = float(st.text_input('Maximum Temperature ('+tempUnit+')', round(simFunctions.tempConversion(tempUnit, 302.1, False), 1)))
                noPoints = st.number_input('Number of Points', 1, None, 4, 1)
                T = numpy.arange(maxTemp, minTemp-(maxTemp-minTemp)/noPoints, -1*(maxTemp-minTemp)/(noPoints-1))
                T = [round(i, 1) for i in T]

                if userGuess == True:
                    c1, c2 = st.columns(2)
                    with c1:
                        minGuessPressure = float(st.text_input("Minimum Guess Pressure ("+pressureUnit+"): ", round(simFunctions.pressureConversion(pressureUnit, 0.995, False),3)))
                    with c2:
                        maxGuessPressure = float(st.text_input("Maximum Guess Pressure ("+pressureUnit+"): ", round(simFunctions.pressureConversion(pressureUnit, 74.291, False),3)))

                    logP = numpy.arange(math.log(maxGuessPressure), math.log(minGuessPressure)-(math.log(maxGuessPressure)-math.log(minGuessPressure))/noPoints, -1*(math.log(maxGuessPressure)-math.log(minGuessPressure))/(noPoints-1))
                    P = [0 for i in range(len(T))]
                    for i in range(len(logP)):
                        T[i] = round(simFunctions.tempConversion(tempUnit, T[i], True), 2)
                        P[i] = math.exp(logP[i])
                else:
                    P = [0 for i in range(len(T))]
                    for i in range(len(T)):
                        T[i] = simFunctions.tempConversion(tempUnit, T[i], True)
                        P[i] = simFunctions.guessPressure(components, moleFractions, T[i])/1E6

            else:
                c1, c2 = st.columns(2)
                with c1:
                    minPressure = float(st.text_input('Minimum Pressure ('+pressureUnit+')', round(simFunctions.pressureConversion(pressureUnit, 0.995, False),3)))
                with c2:
                    maxPressure = float(st.text_input('Maximum Pressure ('+pressureUnit+')', round(simFunctions.pressureConversion(pressureUnit, 74.291, False),3)))
                noPoints = st.number_input('Number of Points', 1, None, 4, 1)
                P = numpy.arange(maxPressure, minPressure-(maxPressure-minPressure)/noPoints, -1*(maxPressure-minPressure)/(noPoints-1))

                if userGuess == True:
                    c1, c2 = st.columns(2)
                    with c1:
                        minGuessTemp = float(st.text_input("Minimum Guess Temperature ("+tempUnit+"): ", round(simFunctions.tempConversion(tempUnit, 243.2, False), 1)))
                    with c2:
                        maxGuessTemp = float(st.text_input("Maximum Guess Pressure ("+tempUnit+"): ", round(simFunctions.tempConversion(tempUnit, 302.1, False), 1)))

                    P = numpy.arange(maxPressure, minPressure-(maxPressure-minPressure)/noPoints, -1*(maxPressure-minPressure)/(noPoints-1))
                    expT = numpy.arange(math.exp(maxGuessTemp/100), math.exp(minGuessTemp/100)-(math.exp(maxGuessTemp/100)-math.exp(minGuessTemp/100))/noPoints, -1*(math.exp(maxGuessTemp/100)-math.exp(minGuessTemp/100))/(noPoints-1))
                    T = [0 for i in range(len(expT))]
                    for i in range(len(expT)):
                        T[i] = round(math.log(expT[i]), 2)*100
                else:
                    T = [0 for i in range(len(P))]
                    for i in range(len(P)):
                        P[i] = simFunctions.pressureConversion(pressureUnit, P[i], True)
                        T[i] = simFunctions.guessTemp(components, moleFractions, P[i]*1E6)

        else:
            if definedVariable == "T":
                T = [round(simFunctions.tempConversion(tempUnit, float(st.text_input('Temperature ('+tempUnit+')', round(simFunctions.tempConversion(tempUnit, 278.1, False), 1))), True), 1)]
                if userGuess == True:
                    P = [round(simFunctions.pressureConversion(pressureUnit, float(st.text_input('Guess Pressure ('+pressureUnit+')', round(simFunctions.pressureConversion(pressureUnit, 4.429, False),3))), True), 3)]
                    if pressureUnit == "MPa":
                        P[0] *= 1E6
                else:
                    P = [simFunctions.guessPressure(components, moleFractions, T[0])]
            else:
                P = [round(simFunctions.pressureConversion(pressureUnit, float(st.text_input('Pressure ('+pressureUnit+')', round(simFunctions.pressureConversion(pressureUnit, 4.429, False), 3))), True), 3)*1E6]
                if userGuess == True:
                    T = [round(simFunctions.tempConversion(tempUnit, float(st.text_input('Guess Temperature ('+tempUnit+')', round(simFunctions.tempConversion(tempUnit, 278.1, False), 1))), True), 1)]
                else:
                    T = [simFunctions.guessTemp(components, moleFractions, P[0])]

    else:
        guessFile = numpy.genfromtxt(csvGuesses, delimiter=',', skip_header=1)
        T = guessFile[:,0]
        P = guessFile[:,1]*1E6
        noPoints = len(T)
        if noPoints > 1:
            calculateRange = True
        else:
            calculateRange = False

    freshWater = st.toggle('Fresh Water', True)
    inhibitorConcs = []
    saltConcs = []
    if freshWater == False:
        if calculateRange == False:
                hydrateType = st.radio("Hydrate Type", ["Pure Methane", "Pure Ethane", "Pure CO2", "Generic Structure I", "Propane", "Generic Structure II"], horizontal=False)
                betaGas = [9.120E-4, 8.165E-4, 9.186E-4, 9.432E-4, 1.058E-3, 8.755E-4][["Pure Methane", "Pure Ethane", "Pure CO2", "Generic Structure I", "Propane", "Generic Structure II"].index(hydrateType)]
        #Inhibitor Concentration input Tables
        salts, inhibitors = simFunctions.getInhibitors()

        c1, c2 = st.columns(2)
        with c1:
            saltConcDf = pd.DataFrame([])
            for i in range(len(salts)):
                saltConcDf = pd.concat([saltConcDf, pd.DataFrame([{'Salt': salts[i], 'Weight Percent': 0.}])], ignore_index=True)
            inputSaltDf = st.data_editor(saltConcDf, hide_index=True)
            saltConcs = inputSaltDf['Weight Percent'].tolist()

        with c2:
            inhibitorConcDf = pd.DataFrame([])
            for i in range(len(inhibitors)):
                inhibitorConcDf = pd.concat([inhibitorConcDf, pd.DataFrame([{'Inhibitor': inhibitors[i], 'Weight Percent': 0.}])], ignore_index=True)
            inputInhibitorDf = st.data_editor(inhibitorConcDf, hide_index=True)
            inhibitorConcs = inputInhibitorDf['Weight Percent'].tolist()

    calculated = False
    if st.button("Calculate"):
        if sum(moleFractions) == 1 and simFunctions.checkMaxConc(inhibitorConcs) == "" and sum(inhibitorConcs)+sum(saltConcs) < 100:
            startTime = time.time()
            eqPressure = [0 for i in range(len(T))]
            if csvGuesses == None:
                if definedVariable == "T":
                    for i in range(len(T)):
                        P[i] = simFunctions.pressureConversion(pressureUnit, P[i], False)
                elif definedVariable == "P":
                    for i in range(len(P)):
                        T[i] = simFunctions.tempConversion(tempUnit, T[i], False)
                        eqTemperature = [0 for i in range(len(P))]

            eqStructure = [0 for i in range(len(T))]
            eqFractions = [0 for i in range(len(T))]
            TInhibited = [0 for i in range(len(T))]
            hydrationNumber = [0 for i in range(len(T))]
            hydrateDensity = [0 for i in range(len(T))]
            if calculateRange == True:
                progressBar = st.progress(0, str(0) + "/" + str(len(T)))
                for i in range(len(T)):
                    if definedVariable == "T":
                        simResult = simFunctions.equilibriumPressure(T[i], P[i], components, moleFractions, saltConcs, inhibitorConcs)
                        eqPressure[i] = simResult[0]
                    elif definedVariable == "P":
                        simResult = simFunctions.equilibriumTemperature(T[i], P[i]*1E6, components, moleFractions, saltConcs, inhibitorConcs)
                        eqTemperature[i] = simResult[0]
                        eqPressure[i] = P[i]
                    eqStructure[i] = simResult[1]
                    eqFractions[i] = [[round(float(simResult[2][0][j]), 4) for j in range(len(simResult[2][0]))], [round(float(simResult[2][1][j]), 4) for j in range(len(simResult[2][1]))]]
                    hydrationNumber[i] = simResult[3]
                    hydrateDensity[i] = simResult[4]
                    with progressBar:
                        st.progress((i+1)/len(T), str(i+1) + "/" + str(len(T)))
            
                if definedVariable == "T":
                    betaGas = simFunctions.betaGas(T, eqPressure)
                elif definedVariable == "P":
                    betaGas = simFunctions.betaGas(eqTemperature, P)
                eqFractions = numpy.array(eqFractions)

                for i in range(len(T)):
                    eqPressure[i] = eqPressure[i]/1E6

                fig, ax = plt.subplots()
                              
                if definedVariable == "T":
                    for i in range(len(T)):
                        eqPressure[i] = simFunctions.pressureConversion(pressureUnit, eqPressure[i], False)
                    "{:.2e}".format(eqPressure[i])
                    P = eqPressure
                elif definedVariable == "P":
                    for i in range(len(P)):
                        eqTemperature[i] = simFunctions.tempConversion(tempUnit, eqTemperature[i], False)
                        eqTemperature[i] = round(eqTemperature[i], 1)
                        T = eqTemperature

                for i in range(len(T)):
                    if freshWater == False:
                        if T[i] >= 273.15:
                            TInhibited[i] = round(simFunctions.HuLeeSum(T[i], saltConcs, inhibitorConcs, betaGas), 1)
                        else:
                            TInhibited[i] = None
                    else:
                        TInhibited[i] = round(T[i], 1)
                    T[i] = simFunctions.tempConversion(tempUnit, T[i], False)
                    if TInhibited[i] != None:
                        TInhibited[i] = round(simFunctions.tempConversion(tempUnit, TInhibited[i], False), 1)

                plt.plot(T, P, '-', label='Fresh Water')

                if freshWater == False:
                    if definedVariable == "T":
                        plt.plot([val for val, condition in zip(TInhibited, TInhibited) if condition is not None], [val for val, condition in zip(eqPressure, TInhibited) if condition is not None], '--', label='Inhibited System')
                    elif definedVariable == "P":
                        plt.plot([val for val, condition in zip(TInhibited, TInhibited) if condition is not None], [val for val, condition in zip(P, TInhibited) if condition is not None], '--', label='Inhibited System')
                    plt.legend(prop={'family': 'Arial'})
                plt.yscale("log")
                plt.xlabel("Temperature ("+tempUnit+")", **{'fontname':'Arial'}, fontsize = 14)
                plt.ylabel("Pressure ("+pressureUnit+")", **{'fontname':'Arial'}, fontsize = 14)
                plt.xticks(**{'fontname':'Arial'}, fontsize = 14)
                plt.yticks(**{'fontname':'Arial'}, fontsize = 14)

                plt.text(0.95, 0.05, "Phases to Flow Research Group", transform=ax.transAxes, fontsize=10, color='gray', alpha=0.8, ha='right', va='bottom', fontweight='bold')
                
                plt.tick_params(axis='both', which='both', direction='in')
                st.pyplot(fig)
                for i in range(len(T)):
                    if definedVariable == "T":
                        if eqPressure[i] < 1 or eqPressure[i] >= 10:
                            eqPressure[i] = f"{eqPressure[i]:.2e}"
                        else:
                            eqPressure[i] = round(eqPressure[i], 2)
                    if definedVariable == "P":
                        eqTemperature[i] = round(eqTemperature[i], 1)
            else:
                if definedVariable == "T":
                    simResult = simFunctions.equilibriumPressure(T[0], P[0], components, moleFractions, saltConcs, inhibitorConcs)
                    eqPressure[0] = simResult[0]/1E6
                    eqTemperature = T
                elif definedVariable == "P":
                    simResult = simFunctions.equilibriumTemperature(T[0], P[0], components, moleFractions, saltConcs, inhibitorConcs)
                    eqTemperature[0] = simResult[0]
                    for i in range(len(P)):
                        P[i] /= 1E6
                    eqPressure = P
                eqStructure[0] = simResult[1]
                eqFractions[0] = [[round(float(simResult[2][0][j]), 4) for j in range(len(simResult[2][0]))], [round(float(simResult[2][1][j]), 4) for j in range(len(simResult[2][1]))]]
                hydrationNumber[0] = simResult[3]
                hydrateDensity[0] = simResult[4]
                eqFractions = numpy.array(eqFractions)
                for i in range(len(T)):
                    if freshWater == False:
                        TInhibited[i] = round(simFunctions.HuLeeSum(eqTemperature[0], saltConcs, inhibitorConcs, betaGas), 1)
                    else:
                        TInhibited[i] = eqTemperature[0]
                    T[i] = round(simFunctions.tempConversion(tempUnit, T[i], True), 1)
                    if TInhibited[i] != None:
                        TInhibited[i] = round(simFunctions.tempConversion(tempUnit, TInhibited[i], True), 1)
                    eqPressure[i] = simFunctions.pressureConversion(pressureUnit, eqPressure[i], True)
                    if eqPressure[i] < 1 or eqPressure[i] >= 10:
                        eqPressure[i] = f"{eqPressure[i]:.2e}"
                    else:
                        eqPressure[i] = round(eqPressure[i], 2)
                    eqTemperature[i] = simFunctions.tempConversion(tempUnit, eqTemperature[i], True)
                    eqTemperature[i] = round(eqTemperature[i], 1)
            endTime = time.time()
            st.text("Time to Complete Calculation: " + str(round(endTime-startTime, 3)) + " seconds")
            if len(T) > 1:
                st.text("per Data Point: " + str(round((endTime-startTime)/noPoints, 3)) + " seconds")
            calculated = True
        else:
            if sum(moleFractions) != 1:
                st.markdown(f":red[Sum of Mole Fractions is Not 1]")
            if simFunctions.checkMaxConc(inhibitorConcs) != "":
                st.markdown(f":red[" + "Inhibitor(s) " + simFunctions.checkMaxConc(inhibitorConcs) + " Exceed(s) Maximum Concentration" + "]")
            if sum(inhibitorConcs)+sum(saltConcs) >= 100:
                st.markdown(f":red[Weight Percent of Inhibitors and Salts Exceeds 100%]")
    if calculated == True:
        if definedVariable == "T":
            P = eqPressure
        elif definedVariable == "P":
            T = eqTemperature
        data = pd.DataFrame({'T ('+tempUnit+')': T, 'Inhibited T ('+tempUnit+')': TInhibited, 'Pressure ('+pressureUnit+')': P, 'Eq. Structure': eqStructure, 'Small Cage Occupancies': eqFractions[:,0].tolist(), 'Large Cage Occupancies': eqFractions[:,1].tolist(), 'Hydration Number': hydrationNumber, 'Hydrate Density (kg/m^3)': hydrateDensity}, [i for i in range(len(T))])
        displayData = pd.DataFrame({'T ('+tempUnit+')': T, 'Inhibited T ('+tempUnit+')': TInhibited, 'Pressure ('+pressureUnit+')': P, 'Eq. Structure': eqStructure}, [i for i in range(len(T))])
        st.dataframe(displayData, hide_index = True)
        st.download_button("Full Data Download", data=data.to_csv(index=False).encode('utf-8'), file_name='data.csv', mime='text/csv')

elif programType == "Minimum Concentration Calculation":
    tempUnit = st.radio("Temperature Unit", ["K", "°C", "°F", "R"], horizontal=True)
    T = float(st.text_input("Fresh Water Equilibrium Temperature ("+tempUnit+"): ", value="280"))
    TDesired = float(st.text_input("Desired Inhibited Equilibrium Temperature ("+tempUnit+"): ", value="275"))
    hydrateType = st.radio("Hydrate Type", ["Pure Methane", "Pure Ethane", "Pure CO2", "Generic Structure I", "Propane", "Generic Structure II"], horizontal=False)
    betaGas = [9.120E-4, 8.165E-4, 9.186E-4, 9.432E-4, 1.058E-3, 8.755E-4][["Pure Methane", "Pure Ethane", "Pure CO2", "Generic Structure I", "Propane", "Generic Structure II"].index(hydrateType)]
    inhibitorType = st.radio("Inhibitor Type", ["Salt", "Liquid Inhibitor"], horizontal=True)
    if inhibitorType == "Liquid Inhibitor":
        salts, inhibitors = simFunctions.getInhibitors()
        inhibitorList = []
        for i in range(len(inhibitors)):
            inhibitorList.append(str(inhibitors[i]))
        inhibitor = st.selectbox('Liquid Inhibitor', inhibitorList)
        inhibitor = inhibitorList.index(inhibitor)
        salt = None
    else:
        inhibitor = "salt"
        salts, inhibitors = simFunctions.getInhibitors()
        saltList = []
        for i in range(len(salts)):
            saltList.append(str(salts[i]))
        salt = st.selectbox('Liquid Inhibitor', saltList)
        salt = saltList.index(salt)
    if st.button("Calculate"):
        T = simFunctions.tempConversion(tempUnit, T, False)
        TDesired = simFunctions.tempConversion(tempUnit, TDesired, False)
        conc = simFunctions.getConcentration(T, TDesired, inhibitor, salt, betaGas, len(inhibitors), len(salts))
        if inhibitor != "salt":
            st.text("Minimum Concentration of " + str(inhibitorList[inhibitor]) + ": " + str(round(conc,1)) + "% w/w")
        else:
            st.text("Minimum Concentration of " + str(saltList[salt]) + ": " + str(round(conc,1)) + "% w/w")

st.header('Credits')
st.markdown('''
            Developed by Karsten Kunneman in collaboration with Prof. Amadeu K. Sum at the Colorado School of Mines \n
            Hydrate model: Klauda-Sandler fugacity model (doi 10.1021/ie000322b) \n
            Inhibhition model: HLS correlation (doi 10.1002/aic.16369)
            \nThis site created with Streamlit''')
st.markdown(f'''<a href="https://github.com/karstenkunneman/Gas-Hydrate-Equilibrium-Calculator">Github Repo</a>''', unsafe_allow_html=True)
