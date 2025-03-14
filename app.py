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

st.title('Gas Hydrate Equilibrium Calculator')
st.caption('Version 2025-03-13')

programType = st.radio("Calculation Type", ["Equilibrium Calculator", "Minimum Concentration Calculator"], horizontal=True)

if programType == "Equilibrium Calculator":
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
    tempUnit = st.radio("Temperature Unit", ["K", "°C", "°F", "R"], horizontal=True)
    pressureUnit = st.radio("Pressure Unit", ["MPa", "bar", "psia", "mmHg"], horizontal=True)

    #If no input file given, take inputs from user through UI
    if csvGuesses == None:
        calculateRange = st.toggle('Calculate Range of Temperatures', False)
        userGuess = st.toggle('Manually Guess Pressures?', False)
        if calculateRange == True:
            c1, c2 = st.columns(2)
            with c1:
                minTemp = float(st.text_input('Minimum Temperature ('+tempUnit+')', 190.15))
            with c2:
                maxTemp = float(st.text_input('Maximum Temperature ('+tempUnit+')', 302.1))
            noPoints = st.number_input('Number of Points', 1, None, 4, 1)
            T = numpy.arange(maxTemp, minTemp-(maxTemp-minTemp)/noPoints, -1*(maxTemp-minTemp)/(noPoints-1))
            T = [round(i, 1) for i in T]

            if userGuess == True:
                c1, c2 = st.columns(2)
                with c1:
                    minGuessPressure = float(st.text_input("Minimum Guess Pressure ("+pressureUnit+"): ", 0.089)) #Pressure in Pa
                with c2:
                    maxGuessPressure = float(st.text_input("Maximum Guess Pressure ("+pressureUnit+"): ", 74.291)) #Pressure in Pa
            else:
                for i in range(len(T)):
                    T[i] = simFunctions.tempConversion(tempUnit, T[i], False)
                minGuessPressure = simFunctions.guessPressure(components, moleFractions, T[len(T)-1])
                maxGuessPressure = simFunctions.guessPressure(components, moleFractions, T[0])

            logP = numpy.arange(math.log(maxGuessPressure), math.log(minGuessPressure)-(math.log(maxGuessPressure)-math.log(minGuessPressure))/noPoints, -1*(math.log(maxGuessPressure)-math.log(minGuessPressure))/(noPoints-1))
            P = [0 for i in range(len(T))]
            for i in range(len(logP)):
                T[i] = round(T[i], 2)
                P[i] = math.exp(logP[i])
        else:
            T = [float(st.text_input('Temperature ('+tempUnit+')', 278.1))]
            if userGuess == True:
                P = [float(st.text_input('Guess Pressure ('+pressureUnit+')', 4.249))]
            else:
                P = [simFunctions.guessPressure(components, moleFractions, T[0])]

    else:
        guessFile = numpy.genfromtxt(csvGuesses, delimiter=',', skip_header=1)
        T = guessFile[:,0]
        P = guessFile[:,1]
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
                saltConcDf = pd.concat([saltConcDf, pd.DataFrame([{'Salt': salts[i], 'Weight Fraction': 0.}])], ignore_index=True)
            inputSaltDf = st.data_editor(saltConcDf, hide_index=True)
            saltConcs = inputSaltDf['Weight Fraction'].tolist()

        with c2:
            inhibitorConcDf = pd.DataFrame([])
            for i in range(len(inhibitors)):
                inhibitorConcDf = pd.concat([inhibitorConcDf, pd.DataFrame([{'Inhibitor': inhibitors[i], 'Weight Fraction': 0.}])], ignore_index=True)
            inputInhibitorDf = st.data_editor(inhibitorConcDf, hide_index=True)
            inhibitorConcs = inputInhibitorDf['Weight Fraction'].tolist()

    calculated = False
    if st.button("Calculate"):
        if sum(moleFractions) == 1 and simFunctions.checkMaxConc(inhibitorConcs) == "" and sum(inhibitorConcs)+sum(saltConcs) < 100:
            startTime = time.time()
            for i in range(len(T)):
                P[i] = simFunctions.pressureConversion(pressureUnit, P[i], False)

            eqPressure = [0 for i in range(len(T))]
            eqStructure = [0 for i in range(len(T))]
            eqFractions = [0 for i in range(len(T))]
            TInhibited = [0 for i in range(len(T))]
            if calculateRange == True:
                progressBar = st.progress(0, str(0) + "/" + str(len(T)))
                for i in range(len(T)):
                    simResult = simFunctions.equilibriumPressure(T[i], P[i], components, moleFractions, saltConcs, inhibitorConcs)
                    eqPressure[i] = simResult[0]
                    eqStructure[i] = simResult[1]
                    eqFractions[i] = [[round(float(simResult[2][0][j]), 4) for j in range(len(simResult[2][0]))], [round(float(simResult[2][1][j]), 4) for j in range(len(simResult[2][1]))]]
                    with progressBar:
                        st.progress((i+1)/len(T), str(i+1) + "/" + str(len(T)))
                betaGas = simFunctions.betaGas(T, eqPressure)
                for i in range(len(T)):
                    if freshWater == False:
                        TInhibited[i] = simFunctions.HuLeeSum(T[i], saltConcs, inhibitorConcs, betaGas)
                    else:
                        TInhibited[i] = T[i]
                eqFractions = numpy.array(eqFractions)
                for i in range(len(T)):
                    T[i] = simFunctions.tempConversion(tempUnit, T[i], True)
                    if TInhibited[i] != None:
                        TInhibited[i] = round(simFunctions.tempConversion(tempUnit, TInhibited[i], True), 1)
                    eqPressure[i] = simFunctions.pressureConversion(pressureUnit, eqPressure[i], True)
                    "{:.2e}".format(eqPressure[i])
                fig, ax = plt.subplots()
                plt.plot(T, eqPressure, '-', label='Fresh Water')
                if freshWater == False:
                    plt.plot([val for val, condition in zip(TInhibited, TInhibited) if condition is not None], [val for val, condition in zip(eqPressure, TInhibited) if condition is not None], '--', label='Inhibited System')
                    plt.legend(prop={'family': 'Arial'})
                plt.yscale("log")
                plt.xlabel("Temperature ("+tempUnit+")", **{'fontname':'Arial'})
                plt.ylabel("Pressure ("+pressureUnit+")", **{'fontname':'Arial'})
                plt.xticks(**{'fontname':'Arial'})
                plt.yticks(**{'fontname':'Arial'})

                plt.text(0.5, 0.5, "Phases to Flow Research Group", transform=ax.transAxes, fontsize=20, color='gray', alpha=0.5, ha='center', va='center', rotation=30)
                
                plt.tick_params(axis='both', which='both', direction='in')
                st.pyplot(fig)
                for i in range(len(T)):
                    if eqPressure[i] < 1 or eqPressure[i] >= 10:
                        eqPressure[i] = f"{eqPressure[i]:.2e}"
                    else:
                        eqPressure[i] = round(eqPressure[i], 2)
                showFracs = st.toggle("Show Fractional Occupancies", False)
            else:
                simResult = simFunctions.equilibriumPressure(T[0], P[0], components, moleFractions, saltConcs, inhibitorConcs)
                eqPressure[0] = simResult[0]
                eqStructure[0] = simResult[1]
                eqFractions[0] = [[round(float(simResult[2][0][j]), 4) for j in range(len(simResult[2][0]))], [round(float(simResult[2][1][j]), 4) for j in range(len(simResult[2][1]))]]
                eqFractions = numpy.array(eqFractions)
                if freshWater == False:
                    TInhibited[i] = round(simFunctions.HuLeeSum(T[0], saltConcs, inhibitorConcs, betaGas), 1)
                else:
                    TInhibited[i] = T[0]
                for i in range(len(T)):
                    T[i] = round(simFunctions.tempConversion(tempUnit, T[i], True), 1)
                    if TInhibited[i] != None:
                        TInhibited[i] = round(simFunctions.tempConversion(tempUnit, TInhibited[i], True), 1)
                    eqPressure[i] = simFunctions.pressureConversion(pressureUnit, eqPressure[i], True)
                    if eqPressure[i] < 1 or eqPressure[i] >= 10:
                        eqPressure[i] = f"{eqPressure[i]:.2e}"
                    else:
                        eqPressure[i] = round(eqPressure[i], 2)
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
        data = pd.DataFrame({'T ('+tempUnit+')': T, 'Inhibited T ('+tempUnit+')': TInhibited, 'Eq. Pressure ('+pressureUnit+')': eqPressure, 'Eq. Structure': eqStructure, 'Small Cage Occupancies': eqFractions[:,0].tolist(), 'Large Cage Occupancies': eqFractions[:,1].tolist()}, [i for i in range(len(T))])
        displayData = pd.DataFrame({'T ('+tempUnit+')': T, 'Inhibited T ('+tempUnit+')': TInhibited, 'Eq. Pressure ('+pressureUnit+')': eqPressure, 'Eq. Structure': eqStructure}, [i for i in range(len(T))])
        st.dataframe(displayData, hide_index = True)
        st.download_button("Full Data Download", data=data.to_csv(index=False).encode('utf-8'), file_name='data.csv', mime='text/csv')

elif programType == "Minimum Concentration Calculator":
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
            Created by Karsten Kunneman and Dr. Amadeu K Sum at the Colorado School of Mines
            \nBased on "Phase behavior of clathrate hydrates: a model for single and multiple 
            gas component hydrates" by Jeffery B. Klauda and Stanley I. Sandler
            \nImplements the Hu-Lee-Sum Correlation for Prediction of Hydrate Phase Equilibria in 
            Mixed Salt and Organic Inhibitor Systems
            \nThis site created with Streamlit''')
st.markdown(f'''<a href="https://github.com/karstenkunneman/Gas-Hydrate-Equilibrium-Calculator">Github Repo</a>''', unsafe_allow_html=True)