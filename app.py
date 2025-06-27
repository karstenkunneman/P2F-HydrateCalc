#This work is licensed under the Creative Commons Attribution-ShareAlike 4.0 International License. To view a copy of this license, visit https://creativecommons.org/licenses/by-sa/4.0/.
#Created by Karsten Kunneman and Amadeu K. Sum at the Colorado School of Mines
#©2025, All Rights Reserved

import simFunctions
import streamlit as st
import pandas as pd
import numpy
import matplotlib.pyplot as plt
from matplotlib import font_manager as fmgr, rcParams
import time
import io

inhibitorData = pd.read_excel('Data.xlsx', sheet_name='Inhibitor Data').to_numpy()

fontPath="static/ARIAL.TTF"
fmgr.fontManager.addfont(fontPath)
fprop = fmgr.FontProperties(fname=fontPath)
fname = fprop.get_name()
rcParams["font.family"] = fname
st.set_page_config(
        page_title="P2F Hydrate Calculator",
        page_icon="thumbnail_P2F_logo(green).png",
    )

@st.cache_data
def equilibriumPressure(temperature, pressure, compounds, moleFractions, saltConcs, inhibitorConcs):
    return simFunctions.equilibriumPressure(temperature, pressure, compounds, moleFractions, saltConcs, inhibitorConcs)

@st.cache_data
def equilibriumTemperature(temperature, pressure, compounds, moleFractions, saltConcs, inhibitorConcs):
    return simFunctions.equilibriumTemperature(temperature, pressure, compounds, moleFractions, saltConcs, inhibitorConcs)

IDs, compounds = simFunctions.getComponents()
componentList = []
for i in range(len(IDs)):
    componentList.append(compounds[i])

c1, c2 = st.columns([2,0.9], gap="small")
with c1:
    st.title('Gas Hydrate Equilibrium Predictions Calculator')
with c2:
    st.image('thumbnail_P2F_logo(green) (FULL).png')

st.caption('Version 1.2.0')

programType = st.radio("Calculation Type", ["Equilibrium Calculation", "Minimum Concentration Calculation"], horizontal=True)

if programType == "Equilibrium Calculation":
    st.caption('NOTE: After selecting "Full Data Download" or "Download Plot", the page will appear to reset. If no changes are made to system parameters, just select "Calculate" again, and you can select other options as desired.')
    #User temperature and guess pressure input file w/ template
    if st.toggle("Upload Temperatures/Pressures", False) == True:
        csvGuesses = st.file_uploader("Upload a .csv file containing temperatures/pressures, (optional) guess pressures/temperatures, and (optional) compositions for each point.", ['csv'])
        csvTemplate = st.download_button("Guess File Template", open("Input Template.csv", encoding='utf-8'), file_name="Input Template.csv")
    else:
        csvGuesses = None
    
    if csvGuesses != None:
        manualComp = st.toggle('Manual Component Input', False)
        userGuess = st.toggle('Guess Pressures from File?', True)

    components = []
    moleFractions = []
    if csvGuesses == None or manualComp == True:
        #Mole Fraction Input Table
        c1, c2 = st.columns(2, gap="medium")
        with c1:
            compDf = pd.DataFrame([])
            compDf = pd.concat([compDf, pd.DataFrame([{'Component': componentList[0], 'Mole Fraction': 1.}])], ignore_index=True)
            for i in range(len(componentList)-1):
                compDf = pd.concat([compDf, pd.DataFrame([{'Component': componentList[i+1], 'Mole Fraction': 0.}])], ignore_index=True)

            inputCompDf = st.data_editor(compDf, hide_index=True)

            moleFracInput = inputCompDf['Mole Fraction'].tolist()

            for i in range(len(componentList)):
                components.append(i + 1)
                moleFractions.append(round(moleFracInput[i],4))

        normalizeFracs = st.toggle("Normalize Mole Fractions", False)

        with c2:
            if normalizeFracs == True:
                nonNormalSum = sum(moleFractions)

                for i in range(len(moleFractions)):
                    moleFractions[i]/=nonNormalSum

                normDf = pd.DataFrame([])
                for i in range(len(componentList)):
                    normDf = pd.concat([normDf, pd.DataFrame([{'Component': componentList[i], 'Normalized Mole Fraction': moleFractions[i]}])], ignore_index=True)

                inputCompDf = st.dataframe(normDf, hide_index=True)

        nonZeroFracs = numpy.nonzero(moleFractions)[0]
        moleFractions = [element for index, element in enumerate(moleFractions) if index in nonZeroFracs]
        components = [element for index, element in enumerate(components) if index in nonZeroFracs]

        st.text("Sum of Mole Fractions: " + str(round(sum(moleFractions),4)))

        if sum(moleFractions) == 1:
            confirmSumFrac = True
        else:
            confirmSumFrac = False

    inhibitorConcs = []
    saltConcs = []
    #Inhibitor Concentration input Tables
    salts, inhibitors = simFunctions.getInhibitors()

    c1, c2 = st.columns(2)
    with c1:
        st.caption("Input Salt Concentration (wt%) based on water amount")
        saltConcDf = pd.DataFrame([])
        for i in range(len(salts)):
            saltConcDf = pd.concat([saltConcDf, pd.DataFrame([{'Salt': salts[i], 'Wt. %': 0.}])], ignore_index=True)
        inputSaltDf = st.data_editor(saltConcDf, hide_index=True)
        saltConcs = inputSaltDf['Wt. %'].tolist()

    with c2:
        st.caption("Organic Inhibitor Concentration (wt%) based on salt aqueous solutions")
        inhibitorConcDf = pd.DataFrame([])
        for i in range(len(inhibitors)):
            inhibitorConcDf = pd.concat([inhibitorConcDf, pd.DataFrame([{'Inhibitor': inhibitors[i], 'Wt. %': 0., 'Max. Conc.': inhibitorData[i][5]}])], ignore_index=True)
        inputInhibitorDf = st.data_editor(inhibitorConcDf, hide_index=True)
        inhibitorConcs = inputInhibitorDf['Wt. %'].tolist()

    #Unit Selector
    tempUnit = st.radio("Temperature Unit", ["K", "°C", "°F"], horizontal=True)
    pressureUnit = st.radio("Pressure Unit", ["MPa", "bar", "psia"], horizontal=True)
    pressureScale = "Standard" #st.radio("Pressure Scale", ["Standard", "Logarithmic"], horizontal=True)

    #Defined Variable Selector
    definedVariable = st.radio("Defined Variable", ["T", "P"], horizontal=True)

    #If no input file given, take inputs from user through UI
    if csvGuesses == None:
        calculateRange = st.toggle('Calculation Range', True)
        userGuess = False
        if calculateRange == True:
            if definedVariable == "T":
                c1, c2 = st.columns(2)
                with c1:
                    minTemp = float(st.text_input('Minimum Temperature ('+tempUnit+')', round(simFunctions.tempConversion(tempUnit, 240.2, False), 1)))
                with c2:
                    maxTemp = float(st.text_input('Maximum Temperature ('+tempUnit+')', round(simFunctions.tempConversion(tempUnit, 302.1, False), 1)))
                noPoints = st.number_input('Number of Points', 1, None, 4, 1)
                T = numpy.arange(minTemp, maxTemp+(maxTemp-minTemp)/noPoints, (maxTemp-minTemp)/(noPoints-1))
                T = [round(i, 1) for i in T]

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


                T = [0 for i in range(len(P))]
                for i in range(len(P)):
                    P[i] = simFunctions.pressureConversion(pressureUnit, P[i], True)
                    T[i] = simFunctions.guessTemp(components, moleFractions, P[i]*1E6)

        else:
            hydrateType = st.radio("Hydrate Type", ["Pure Methane", "Pure Ethane", "Pure CO2", "Generic Structure I", "Propane", "Generic Structure II"], horizontal=False)
            betaGas = [9.120E-4, 8.165E-4, 9.186E-4, 9.432E-4, 1.058E-3, 8.755E-4][["Pure Methane", "Pure Ethane", "Pure CO2", "Generic Structure I", "Propane", "Generic Structure II"].index(hydrateType)]
            if definedVariable == "T":
                T = [round(simFunctions.tempConversion(tempUnit, float(st.text_input('Temperature ('+tempUnit+')', round(simFunctions.tempConversion(tempUnit, 278.1, False), 1))), True), 1)]
                P = [simFunctions.guessPressure(components, moleFractions, T[0])]
            else:
                P = [round(simFunctions.pressureConversion(pressureUnit, float(st.text_input('Pressure ('+pressureUnit+')', round(simFunctions.pressureConversion(pressureUnit, 4.429, False), 3))), True), 3)*1E6]
                T = [simFunctions.guessTemp(components, moleFractions, P[0])]

    else:
        guessFile = numpy.genfromtxt(csvGuesses, delimiter=',', skip_header=1)
        T = guessFile[:,0]
        P = guessFile[:,1]
        noPoints = len(T)
        if noPoints > 1:
            calculateRange = True
        else:
            calculateRange = False
        if manualComp == False:
            for i in range(len(T)):
                pointComponents = []
                pointMoleFractions = []
                for j in range(len(IDs)):
                    if guessFile[:,2+j][i] != 0 and numpy.isnan(guessFile[:,2+j][i]) == False:
                        pointComponents.append(j+1)
                        pointMoleFractions.append(guessFile[:,2+j][i])
                components += [pointComponents]
                moleFractions += [pointMoleFractions]
            for i in range(len(components)):
                if round(sum(moleFractions[i]),4) == 1 or confirmSumFrac == True:
                    confirmSumFrac = True
                else:
                    st.markdown(f":red[Sum of component fractions in row " + str(i+1)+ " do not equal 1")
                    confirmSumFrac = False
        else:
            components = [components for i in range(noPoints)]
            moleFractions = [moleFractions for i in range(noPoints)]
        if userGuess == False:
            P = [0 for i in range(noPoints)]
            for i in range(noPoints):
                T[i] = simFunctions.tempConversion(tempUnit, T[i], True)
                P[i] = simFunctions.guessPressure(components[i], moleFractions[i], T[i])/1E6

    calculated = False
    if st.button("Calculate"):
        if confirmSumFrac == True and simFunctions.checkMaxConc(inhibitorConcs) == "" and sum(inhibitorConcs)+sum(saltConcs) < 100:
            startTime = time.time()
            eqPressure = [0 for i in range(len(T))]
            if csvGuesses == None:
                if definedVariable == "T":
                    for i in range(len(T)):
                        if userGuess == True:
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
            eqPhase = [0 for i in range(len(T))]
            if calculateRange == True:
                progressBar = st.progress(0, str(0) + "/" + str(len(T)))
                for i in range(len(T)):
                    if definedVariable == "T":
                        try:
                            if manualComp == True and csvGuesses == None:
                                simResult = equilibriumPressure(T[i], P[i]*1E6, components, moleFractions, saltConcs, inhibitorConcs)
                            else:
                                simResult = equilibriumPressure(T[i], P[i]*1E6, components[i], moleFractions[i], saltConcs, inhibitorConcs)
                        except:
                            simResult = equilibriumPressure(T[i], P[i]*1E6, components, moleFractions, saltConcs, inhibitorConcs)
                        eqPressure[i] = simResult[0]
                    elif definedVariable == "P":
                        try:
                            if manualComp == True and csvGuesses == None:
                                simResult = equilibriumTemperature(T[i], P[i]*1E6, components, moleFractions, saltConcs, inhibitorConcs)
                            else:
                                simResult = equilibriumTemperature(T[i], P[i]*1E6, components[i], moleFractions[i], saltConcs, inhibitorConcs)
                        except:
                            simResult = equilibriumTemperature(T[i], P[i]*1E6, components, moleFractions, saltConcs, inhibitorConcs)
                        eqTemperature[i] = simResult[0]
                        eqPressure[i] = P[i]
                    eqStructure[i] = simResult[1]
                    eqFractions[i] = [[round(float(simResult[2][0][j]), 4) for j in range(len(simResult[2][0]))], [round(float(simResult[2][1][j]), 4) for j in range(len(simResult[2][1]))]]
                    hydrationNumber[i] = simResult[3]
                    hydrateDensity[i] = simResult[4]
                    eqPhase[i] = simResult[5]
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
                    P = eqPressure
                elif definedVariable == "P":
                    for i in range(len(P)):
                        eqTemperature[i] = simFunctions.tempConversion(tempUnit, eqTemperature[i], False)
                        eqTemperature[i] = round(eqTemperature[i], 1)
                        T = eqTemperature

                for i in range(len(T)):
                    if sum(inhibitorConcs) + sum(saltConcs) != 0:
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

                if sum(inhibitorConcs) + sum(saltConcs) != 0:
                    if definedVariable == "T":
                        plt.plot([val for val, condition in zip(TInhibited, TInhibited) if condition is not None], [val for val, condition in zip(eqPressure, TInhibited) if condition is not None], '--', label='Inhibited System')
                    elif definedVariable == "P":
                        plt.plot([val for val, condition in zip(TInhibited, TInhibited) if condition is not None], [val for val, condition in zip(P, TInhibited) if condition is not None], '--', label='Inhibited System')
                    plt.legend()
                if pressureScale == "Logarithmic":
                    plt.yscale("log")
                plt.xlabel("Temperature ("+tempUnit+")", fontsize = 14)
                plt.ylabel("Pressure ("+pressureUnit+")", fontsize = 14)
                plt.xticks(fontsize = 14)
                plt.yticks(fontsize = 14)

                plt.text(0.95, 0.05, "P2F Lab", transform=ax.transAxes, fontsize=10, color='gray', alpha=0.8, ha='right', va='bottom', fontweight='bold')
                
                plt.tick_params(axis='both', which='both', direction='in')
                st.pyplot(fig)
                for i in range(len(T)):
                    if definedVariable == "T":
                        eqPressure[i] = round(eqPressure[i], 2)
                    if definedVariable == "P":
                        eqTemperature[i] = round(eqTemperature[i], 1)
            else:
                if definedVariable == "T":
                    simResult = equilibriumPressure(T[0], P[0], components, moleFractions, saltConcs, inhibitorConcs)
                    eqPressure[0] = simResult[0]/1E6
                    eqTemperature = T
                elif definedVariable == "P":
                    simResult = equilibriumTemperature(T[0], P[0], components, moleFractions, saltConcs, inhibitorConcs)
                    eqTemperature[0] = simResult[0]
                    for i in range(len(P)):
                        P[i] /= 1E6
                    eqPressure = P
                eqStructure[0] = simResult[1]
                eqFractions[0] = [[round(float(simResult[2][0][j]), 4) for j in range(len(simResult[2][0]))], [round(float(simResult[2][1][j]), 4) for j in range(len(simResult[2][1]))]]
                hydrationNumber[0] = simResult[3]
                hydrateDensity[0] = simResult[4]
                eqPhase[0] = simResult[5]
                eqFractions = numpy.array(eqFractions)
                for i in range(len(T)):
                    if sum(inhibitorConcs) + sum(saltConcs) != 0:
                        TInhibited[i] = round(simFunctions.HuLeeSum(eqTemperature[0], saltConcs, inhibitorConcs, betaGas), 1)
                    else:
                        TInhibited[i] = eqTemperature[0]
                    T[i] = round(simFunctions.tempConversion(tempUnit, T[i], True), 1)
                    if TInhibited[i] != None:
                        TInhibited[i] = round(simFunctions.tempConversion(tempUnit, TInhibited[i], True), 1)
                    eqPressure[i] = simFunctions.pressureConversion(pressureUnit, eqPressure[i], True)
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
        if sum(inhibitorConcs) + sum(saltConcs) != 0:
            displayData = pd.DataFrame({'Temperature ('+tempUnit+')': T, 'Inhibited T ('+tempUnit+')': TInhibited, 'Pressure ('+pressureUnit+')': P, 'Eq. Structure': eqStructure}, [i for i in range(len(T))])
            data = pd.DataFrame({'T ('+tempUnit+')': T, 'Inhibited T ('+tempUnit+')': TInhibited, 'Pressure ('+pressureUnit+')': P, 'Eq. Structure': eqStructure, 'Small Cage Occupancies': eqFractions[:,0].tolist(), 'Large Cage Occupancies': eqFractions[:,1].tolist(), 'Hydration Number': hydrationNumber, 'Hydrate Density (kg/m^3)': hydrateDensity, 'Phase Equilibrium': eqPhase}, [i for i in range(len(T))])
        else:
            displayData = pd.DataFrame({'Temperature ('+tempUnit+')': T, 'Pressure ('+pressureUnit+')': P, 'Eq. Structure': eqStructure}, [i for i in range(len(T))])
            data = pd.DataFrame({'T ('+tempUnit+')': T, 'Pressure ('+pressureUnit+')': P, 'Eq. Structure': eqStructure, 'Small Cage Occupancies': eqFractions[:,0].tolist(), 'Large Cage Occupancies': eqFractions[:,1].tolist(), 'Hydration Number': hydrationNumber, 'Hydrate Density (kg/m^3)': hydrateDensity, 'Phase Equilibrium': eqPhase}, [i for i in range(len(T))])
        st.dataframe(displayData, hide_index = True)
        c1, c2 = st.columns(2)
        with c1:
            st.download_button("Save Full Dataset", data=data.to_csv(index=False).encode('utf-8'), file_name='data.csv', mime='text/csv')
        with c2:
            img = io.BytesIO()
            plt.savefig(img, format='png')
            downloadButton = st.download_button(label="Save Plot", data=img, file_name="Plot.png", mime="image/png")
        st.text("Full dataset includes temperature, inhibited temperature, pressure, structure, small and large cage occupancy, hydration number, hydrate density, and phase for each point")

        if st.button("Reset"):
            st.resetpage()

elif programType == "Minimum Concentration Calculation":
    tempUnit = st.radio("Temperature Unit", ["K", "°C", "°F"], horizontal=True)
    T = float(st.text_input("Fresh Water Equilibrium Temperature ("+tempUnit+"): ", value="280"))
    TDesired = float(st.text_input("Desired Inhibited Equilibrium Temperature ("+tempUnit+"): ", value="275"))

    components = []
    moleFractions = []

    #Mole Fraction Input Table
    c1, c2 = st.columns(2)
    with c1:
        compDf = pd.DataFrame([])
        compDf = pd.concat([compDf, pd.DataFrame([{'Component': componentList[0], 'Mole Fraction': 1.}])], ignore_index=True)
        for i in range(len(componentList)-1):
            compDf = pd.concat([compDf, pd.DataFrame([{'Component': componentList[i+1], 'Mole Fraction': 0.}])], ignore_index=True)

        inputCompDf = st.data_editor(compDf, hide_index=True)

        moleFracInput = inputCompDf['Mole Fraction'].tolist()

        for i in range(len(componentList)):
            components.append(i + 1)
            moleFractions.append(round(moleFracInput[i],4))

    normalizeFracs = st.toggle("Normalize Mole Fractions", False)

    with c2:
        if normalizeFracs == True:
            nonNormalSum = sum(moleFractions)

            for i in range(len(moleFractions)):
                moleFractions[i]/=nonNormalSum

            normDf = pd.DataFrame([])
            for i in range(len(componentList)):
                normDf = pd.concat([normDf, pd.DataFrame([{'Component': componentList[i], 'Normalized Mole Fraction': moleFractions[i]}])], ignore_index=True)

            inputCompDf = st.dataframe(normDf, hide_index=True)

    nonZeroFracs = numpy.nonzero(moleFractions)[0]
    moleFractions = [element for index, element in enumerate(moleFractions) if index in nonZeroFracs]
    components = [element for index, element in enumerate(components) if index in nonZeroFracs]

    st.text("Sum of Mole Fractions: " + str(round(sum(moleFractions),4)))

    if sum(moleFractions) == 1:
        confirmSumFrac = True
    else:
        confirmSumFrac = False

    #hydrateType = st.radio("Hydrate Type", ["Pure Methane", "Pure Ethane", "Pure CO2", "Generic Structure I", "Propane", "Generic Structure II"], horizontal=False)
    #betaGas = [9.120E-4, 8.165E-4, 9.186E-4, 9.432E-4, 1.058E-3, 8.755E-4][["Pure Methane", "Pure Ethane", "Pure CO2", "Generic Structure I", "Propane", "Generic Structure II"].index(hydrateType)]
    inhibitorType = st.radio("Inhibitor Type", ["Salt", "Organic Inhibitor"], horizontal=True)
    if inhibitorType == "Organic Inhibitor":
        salts, inhibitors = simFunctions.getInhibitors()
        inhibitorList = []
        for i in range(len(inhibitors)):
            inhibitorList.append(str(inhibitors[i]))
        inhibitor = st.selectbox('Organic Inhibitor', inhibitorList)
        inhibitor = inhibitorList.index(inhibitor)
        salt = None
    else:
        inhibitor = "salt"
        salts, inhibitors = simFunctions.getInhibitors()
        saltList = []
        for i in range(len(salts)):
            saltList.append(str(salts[i]))
        salt = st.selectbox('Salt Type', saltList)
        salt = saltList.index(salt)

    if st.button("Calculate"):
        if confirmSumFrac != True:
            st.markdown(f":red[Sum of component fractions in row " + str(i+1)+ " do not equal 1")
            st.reset()
        T = simFunctions.tempConversion(tempUnit, T, False)
        Tlist = [T-0.5, T, T+0.5]

        guessP = [simFunctions.guessPressure(components, moleFractions, Tlist[0]), simFunctions.guessPressure(components, moleFractions, Tlist[1]), simFunctions.guessPressure(components, moleFractions, Tlist[2])]

        P = [equilibriumPressure(Tlist[0], guessP[0], components, moleFractions, [], [])[0],
             equilibriumPressure(Tlist[1], guessP[1], components, moleFractions, [], [])[0],
             equilibriumPressure(Tlist[2], guessP[2], components, moleFractions, [], [])[0]]
        
        betaGas = simFunctions.betaGas(Tlist, P)

        TDesired = simFunctions.tempConversion(tempUnit, TDesired, False)
        if inhibitor != "salt":
            conc = simFunctions.getConcentration(T, TDesired, inhibitor, salt, betaGas, 1, 0)
            st.text("Minimum Concentration of " + str(inhibitorList[inhibitor]) + ": " + str(round(conc,1)) + "% w/w")
        else:
            conc = simFunctions.getConcentration(T, TDesired, inhibitor, salt, betaGas, 0, 1)
            st.text("Minimum Concentration of " + str(saltList[salt]) + ": " + str(round(conc,1)) + "% w/w")

st.caption("Disclaimer: The model and predictions have been tested and verified to be accurate based on extensive comparison with available literature data. However, this web app is provided ""as is"" and ""as available"" without any warranties of any kind, either express or implied, including, but not limited to, implied warranties of use, merchantability, fitness for a particular purpose, and non-infringement.")

st.header('Credits')
st.markdown('''
            Developed by Karsten Kunneman in collaboration with Prof. Amadeu K. Sum at the Colorado School of Mines
            \nHydrate model: Klauda-Sandler fugacity model [[doi:10.1021/ie000322b]](https://doi.org/10.1021/ie000322b)
            \nInhibition model: HLS correlation [[doi:10.1002/aic.16369]](https://doi.org/10.1002/aic.16369)
            \nThis site is created with Streamlit''')
st.markdown(f'''<a href="https://github.com/karstenkunneman/Gas-Hydrate-Equilibrium-Calculator">Github Repo</a>''', unsafe_allow_html=True)
