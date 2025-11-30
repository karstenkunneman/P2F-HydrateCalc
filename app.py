#This work is licensed under the Creative Commons Attribution-ShareAlike 4.0 International License. To view a copy of this license, visit https://creativecommons.org/licenses/by-sa/4.0/.
#Created by Karsten Kunneman and Amadeu K. Sum at the Colorado School of Mines
#©2025, All Rights Reserved

import simFunctions
import p2f_HydrateCalcLib.core as core
import p2f_HydrateCalcLib.model as model
import streamlit as st
import pandas as pd
import numpy
import matplotlib.pyplot as plt
from matplotlib import font_manager as fmgr, rcParams
import time
import io
import base64
import math

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
def equilibrium(definedVariable, components, moleFractions, T=None, P=None):
    simulation = model.KlaudaSandler2003(components, moleFractions, definedVariable, T, P)
    eqTemperature = simulation.temperature
    eqPressure = simulation.pressure
    eqStructure = simulation.eqStructure
    eqFractions = simulation.eqFrac
    hydrationNumber = simulation.hydrationNumber
    hydrateDensity = [simulation.hydrateDensity, simulation.storageDensity]
    eqPhase = simulation.eqPhase
    freezingPoint = simulation.freezingPoint
    return eqTemperature, eqPressure, eqStructure, eqFractions, hydrationNumber, hydrateDensity, eqPhase, freezingPoint

IDs, compounds = simFunctions.getComponents()
componentList = []
for i in range(len(IDs)):
    componentList.append(compounds[i])

st.image("Mines-horiz-white.png",width=400)

c1, c2 = st.columns([2,0.9], gap="small")
with c1:
    st.title('Gas Hydrate Equilibrium Predictions Calculator')
with c2:
    #st.markdown("[![Phases to Flow Laboratory](thumbnail_P2F_logo(green) (FULL).png)](people.mines.edu/asum/)", unsafe_allow_html=True)
    st.markdown(
    """<a href="https://people.mines.edu/asum/">
    <img src="data:image/png;base64,{}" width="462">
    </a>""".format(
        base64.b64encode(open("thumbnail_P2F_logo(green) (FULL).png", "rb").read()).decode()
    ), unsafe_allow_html=True)

st.caption('Version 1.4.0')

programType = st.radio("Calculation Type", ["Equilibrium Calculation", "Minimum Concentration Calculation"], horizontal=True)

if programType == "Equilibrium Calculation":
    def reset_widgets():
        for key in st.session_state.items():
            del st.session_state[key]

    #User temperature and guess pressure input file w/ template
    if st.toggle("Upload Temperatures/Pressures", False) == True:
        csvGuesses = st.file_uploader("Upload a .csv file containing temperatures/pressures, (optional) guess pressures/temperatures, and (optional) compositions for each point.", ['csv'])
        csvTemplate = st.download_button("Guess File Template", open("Input Template.csv", encoding='utf-8'), file_name="Input Template.csv")
        manualComp = False
    else:
        csvGuesses = None
        manualComp = True
    
    if csvGuesses != None:
        manualComp = st.toggle('Manual Component Input', False)
        userGuess = st.toggle('Guess from File?', False)

    components = []
    moleFractions = []
    if csvGuesses == None or manualComp == True:
        #Mole Fraction Input Table
        massFraction = st.toggle("Input Mass Fractions", False)
        c1, c2, c3 = st.columns([1, 0.9, 1.1], gap="small")
        with c1:
            if massFraction == False:
                compDf = pd.DataFrame([])
                compDf = pd.concat([compDf, pd.DataFrame([{'Component': componentList[0], 'Mole Fraction': 1.}])], ignore_index=True)
                for i in range(len(componentList)-1):
                    compDf = pd.concat([compDf, pd.DataFrame([{'Component': componentList[i+1], 'Mole Fraction': 0.}])], ignore_index=True)

                inputCompDf = st.data_editor(compDf, hide_index=True, column_config={
                    "Component": st.column_config.TextColumn("Component", disabled=True),
                    "Mole Fraction": st.column_config.NumberColumn("Mole Fraction"),
                })

                moleFracInput = inputCompDf['Mole Fraction'].tolist()

                for i in range(len(componentList)):
                    components.append(i + 1)
                    moleFractions.append(round(moleFracInput[i],4))

                moleFractions = numpy.nan_to_num(moleFractions, nan=0)
                    
            else:
                massFractions = []
                compDf = pd.DataFrame([])
                compDf = pd.concat([compDf, pd.DataFrame([{'Component': componentList[0], 'Mass Fraction': 1.}])], ignore_index=True)
                for i in range(len(componentList)-1):
                    compDf = pd.concat([compDf, pd.DataFrame([{'Component': componentList[i+1], 'Mass Fraction': 0.}])], ignore_index=True)

                inputCompDf = st.data_editor(compDf, hide_index=True, column_config={
                    "Component": st.column_config.TextColumn("Component", disabled=True),
                    "Mass Fraction": st.column_config.NumberColumn("Mass Fraction"),
                })

                massFracInput = inputCompDf['Mass Fraction'].tolist()

                for i in range(len(componentList)):
                    components.append(i + 1)
                    massFractions.append(round(massFracInput[i],4))

        c21, c22 = st.columns(2)

        with c21:
            normalizeFracs = st.toggle("Normalize Mole Fractions", False)

        with c2:
            if massFraction == True:
                moleFractions = simFunctions.massToMolFrac(compounds, massFractions)

                normDf = pd.DataFrame([])
                for i in range(len(componentList)):
                    normDf = pd.concat([normDf, pd.DataFrame([{'Component': componentList[i], 'Mole Fraction': moleFractions[i]}])], ignore_index=True)

                inputCompDf = st.dataframe(normDf, hide_index=True)

        with c3:
            if normalizeFracs == True:
                nonNormalSum = sum(moleFractions)

                for i in range(len(moleFractions)):
                    try:
                        moleFractions[i]/=nonNormalSum
                    except:
                        pass

                normDf = pd.DataFrame([])
                for i in range(len(componentList)):
                    normDf = pd.concat([normDf, pd.DataFrame([{'Component': componentList[i], 'Normalized Mole Frac.': moleFractions[i]}])], ignore_index=True)

                inputCompDf = st.dataframe(normDf, hide_index=True)

        nonZeroFracs = numpy.nonzero(numpy.nan_to_num(moleFractions, nan=0))[0]
        moleFractions = [element for index, element in enumerate(moleFractions) if index in nonZeroFracs]
        components = [element for index, element in enumerate(components) if index in nonZeroFracs]

        with c22:
            st.text("Sum of Mole Fractions: " + str(round(sum(moleFractions),4)))

        if sum(moleFractions) == 1:
            confirmSumFrac = True
        else:
            confirmSumFrac = False

    inhibitorConcs = []
    saltConcs = []
    #Inhibitor Concentration input Tables
    salts, inhibitors = simFunctions.getInhibitors()

    saltConcs = numpy.zeros(len(salts))
    inhibitorConcs = numpy.zeros(len(inhibitors))

    notFreshWater = st.toggle("Add Thermodynamic Inhbitors?", False)

    if notFreshWater:
        c1, c2 = st.columns(2, gap="medium")
        with c1:
            st.caption("Input Salt Concentration (wt%) based on water amount")
            saltConcDf = pd.DataFrame([])
            for i in range(len(salts)):
                saltConcDf = pd.concat([saltConcDf, pd.DataFrame([{'Salt': salts[i], 'Wt. %': 0.}])], ignore_index=True)
            inputSaltDf = st.data_editor(saltConcDf, hide_index=True, column_config={
                "Salt": st.column_config.TextColumn("Salt", disabled=True),
                "Wt. %": st.column_config.NumberColumn("Wt. %"),
            })
            saltConcs = inputSaltDf['Wt. %'].tolist()

        with c2:
            st.caption("Organic Inhibitor Concentration (wt%) based on salt aqueous solutions")
            inhibitorConcDf = pd.DataFrame([])
            for i in range(len(inhibitors)):
                inhibitorConcDf = pd.concat([inhibitorConcDf, pd.DataFrame([{'Inhibitor': inhibitors[i], 'Wt. %': 0., 'Max. Conc.': core.inhibitorData[i][5]}])], ignore_index=True)
            inputInhibitorDf = st.data_editor(inhibitorConcDf, hide_index=True, column_config={
                "Inhibitor": st.column_config.TextColumn("Inhibitor", disabled=True),
                "Wt. %": st.column_config.NumberColumn("Wt. %"),
                "Max. Conc.": st.column_config.TextColumn("Max. Conc.", disabled=True)
            })
        inhibitorConcs = inputInhibitorDf['Wt. %'].tolist()

    c1, c2, c3 = st.columns(3)
    with c1:
        tempUnit = st.radio("Temperature Unit", ["K", "°C", "°F"], horizontal=True)
    with c2:
        pressureUnit = st.radio("Pressure Unit", ["MPa", "bar", "psia"], horizontal=True, index=0)
    with c3:
        definedVariable = st.radio("Defined Variable", ["T", "P"], horizontal=True)
    pressureScale = "Standard" #st.radio("Pressure Scale", ["Standard", "Logarithmic"], horizontal=True)

    #If no input file given, take inputs from user through UI
    if csvGuesses == None:
        calculateRange = st.toggle('Calculation Range', True)
        userGuess = False
        if calculateRange == True:
            if definedVariable == "T":
                c1, c2 = st.columns(2)
                with c1:
                    minTemp = float(st.text_input('Minimum Temperature ('+tempUnit+')', round(simFunctions.tempConversion(tempUnit, 262.4, False), 1)))
                with c2:
                    maxTemp = float(st.text_input('Maximum Temperature ('+tempUnit+')', round(simFunctions.tempConversion(tempUnit, 283.8, False), 1)))
                noPoints = st.number_input('Number of Points', 1, None, 8, 1)
                T = numpy.arange(minTemp, maxTemp+(maxTemp-minTemp)/noPoints, (maxTemp-minTemp)/(noPoints-1))

                components = [components for i in range(len(T))]
                moleFractions = [moleFractions for i in range(len(T))]

                P = [0 for i in range(len(T))]
                for i in range(len(T)):
                    T[i] = simFunctions.tempConversion(tempUnit, T[i], True)
                    P[i] = core.guessPressure(core.getComponentData(components[i]), moleFractions[i], T[i])/1E6

            else:
                c1, c2 = st.columns(2)
                with c1:
                    minPressure = float(st.text_input('Minimum Pressure ('+pressureUnit+')', round(simFunctions.pressureConversion(pressureUnit, 1.798, False),3)))
                with c2:
                    maxPressure = float(st.text_input('Maximum Pressure ('+pressureUnit+')', round(simFunctions.pressureConversion(pressureUnit, 7.41, False),3)))
                noPoints = st.number_input('Number of Points', 1, None, 8, 1)
                P = numpy.arange(maxPressure, minPressure-(maxPressure-minPressure)/noPoints, -1*(maxPressure-minPressure)/(noPoints-1))

                components = [components for i in range(len(P))]
                moleFractions = [moleFractions for i in range(len(P))]

                T = [None for i in range(len(P))]
                for i in range(len(P)):
                    P[i] = simFunctions.pressureConversion(pressureUnit, P[i], True)

        else:
            noPoints = 1

            T = simFunctions.tempConversion(tempUnit, float(st.text_input('Temperature ('+tempUnit+')', round(simFunctions.tempConversion(tempUnit, 278.1, False),1))), True)

            if confirmSumFrac != True:
                st.markdown(f":red[Sum of component fractions in row " + str(i+1)+ " do not equal 1")
                st.reset()
            T = simFunctions.tempConversion(tempUnit, T, False)
            T = [T-0.5, T, T+0.5]
            P = [None, None, None]         

    else:
        guessFile = numpy.genfromtxt(csvGuesses, delimiter=',', skip_header=1)
        if definedVariable == "T":
            T = guessFile[:,0]
            noPoints = len(T)
        else:
            P = guessFile[:,1]
            noPoints = len(P)
        if noPoints > 1:
            calculateRange = True
        else:
            calculateRange = False
        if manualComp == False:
            for i in range(noPoints):
                pointComponents = []
                pointMoleFractions = []
                for j in range(len(IDs)):
                    if guessFile[:,2+j][i] != 0 and numpy.isnan(guessFile[:,2+j][i]) == False:
                        pointComponents.append(j+1)
                        pointMoleFractions.append(guessFile[:,2+j][i])
                components += [pointComponents]
                moleFractions += [pointMoleFractions]
                if any(item is not None for item in components[i]):
                    pass
                else:
                    components[i] = components[i-1]
                if any(item is not None for item in moleFractions[i]):
                    pass
                else:
                    moleFractions[i] = moleFractions[i-1]
            for i in range(len(components)):
                if round(sum(moleFractions[i]),4) == 1:
                    confirmSumFrac = True
                else:
                    st.markdown(f":red[Sum of component fractions in row " + str(i+1)+ " do not equal 1")
                    confirmSumFrac = False
        else:
            components = [components for i in range(noPoints)]
            moleFractions = [moleFractions for i in range(noPoints)]

        if definedVariable == "T":
            P = guessFile[:,1]
            if userGuess != True:
                P = [0 for i in range(noPoints)]
                for i in range(noPoints):
                    P[i] = core.guessPressure(core.getComponentData(components[i]), moleFractions[i], T[i])/1E6
            else:
                for i in range(noPoints):
                    if P[i] is None or math.isnan(P[i]):
                        st.markdown(f":red[Guess pressure(s) missing, turn off manual guessing if this is intentional]")
                        break
                    if T[i] is None or math.isnan(T[i]):
                        st.markdown(f":red[Temperature(s) missing, change defined variable if this is intentional]")
                        break
                
        else:
            T = guessFile[:,0]
            if userGuess != True:
                T = [0 for i in range(noPoints)]
                for i in range(noPoints):
                    T[i] = simFunctions.guessTemp(components[i], moleFractions[i], P[i]*1E6)
            else:
                for i in range(noPoints):
                    if T[i] is None or math.isnan(T[i]):
                        st.markdown(f":red[Guess temperature(s) missing, turn off manual guessing if this is intentional]")
                        break
                    if P[i] is None or math.isnan(P[i]):
                        st.markdown(f":red[Pressure(s) missing, change defined variable if this is intentional]")
                        break

    calculated = False
    col1,col2,col3 = st.columns(3)
    with col1:
        TDecimals = st.number_input("Decimal Points in Temperature", value=1, min_value = 0)
    with col2:
        PDecimals = st.number_input("Decimal Points in Pressure", value=2, min_value = 0)

    if st.button("Calculate"):
        if confirmSumFrac == True and simFunctions.checkMaxConc(inhibitorConcs) == "" and sum(inhibitorConcs)+sum(saltConcs) < 100:
            startTime = time.time()
            eqPressure = [0 for i in range(len(T))]
            eqTemperature = [0 for i in range(len(P))]

            simResult = [[] for i in range(len(T))]
            eqStructure = [0 for i in range(len(T))]
            eqFractions = [0 for i in range(len(T))]
            TInhibited = [0 for i in range(len(T))]
            hydrationNumber = [0 for i in range(len(T))]
            hydrateDensity = [0 for i in range(len(T))]
            eqPhase = [0 for i in range(len(T))]
            if calculateRange == True:
                progressBar = st.progress(0, str(0) + "/" + str(len(T)))
                for i in range(len(T)):
                    simResult[i] =  equilibrium(definedVariable, components[i], moleFractions[i], T[i], P[i]*1E6)
                    eqTemperature[i] = simResult[i][0]
                    eqPressure[i] = simResult[i][1]
                    eqStructure[i] = simResult[i][2]
                    eqFractions[i] = [[round(float(simResult[i][3][0][j]), 4) for j in range(len(simResult[i][3][0]))], [round(float(simResult[i][3][1][j]), 4) for j in range(len(simResult[i][3][1]))]]
                    hydrationNumber[i] = simResult[i][4]
                    hydrateDensity[i] = simResult[i][5]
                    eqPhase[i] = simResult[i][6]
                    freezingPoint = simResult[i][7]
                    with progressBar:
                        st.progress((i+1)/len(T), str(i+1) + "/" + str(len(T)))
            
                betaGas = core.betaGas(eqTemperature, eqPressure)
                eqFractions = numpy.array(eqFractions)

                for i in range(len(T)):
                    eqPressure[i] = eqPressure[i]/1E6

                fig, ax = plt.subplots()
                              
                P = eqPressure
                T = eqTemperature

                for i in range(len(T)):
                    if sum(inhibitorConcs) + sum(saltConcs) != 0:
                        if T[i] >= 273.15:
                            TInhibited[i] = round(core.HuLeeSum(T[i], saltConcs, inhibitorConcs, betaGas, freezingPoint), TDecimals)
                        else:
                            TInhibited[i] = None
                    else:
                        TInhibited[i] = round(T[i], TDecimals)
                    T[i] = simFunctions.tempConversion(tempUnit, T[i], False)
                    P[i] = simFunctions.pressureConversion(pressureUnit, P[i], False)
                    if TInhibited[i] != None:
                        TInhibited[i] = round(simFunctions.tempConversion(tempUnit, TInhibited[i], False), TDecimals)

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
                        eqPressure[i] = round(eqPressure[i], PDecimals)
                    if definedVariable == "P":
                        eqTemperature[i] = round(eqTemperature[i], TDecimals)
            else:
                if definedVariable == "T":
                    simResults = [equilibrium("T", components, moleFractions, T=T[0], P=None),
                        equilibrium("T", components, moleFractions, T=T[1], P=None),
                        equilibrium("T", components, moleFractions, T=T[2], P=None)]
                else:
                    P = simFunctions.pressureConversion(pressureUnit, float(st.text_input('Pressure ('+pressureUnit+')', round(simFunctions.pressureConversion(pressureUnit, 4.429, False),1))), True)*1E6
                    Plist = [P*.99, P, P*1.01]

                    simResults = [equilibrium("P", components, moleFractions, T=None, P=Plist[0]),
                        equilibrium("P", components, moleFractions, T=None, P=Plist[1]),
                        equilibrium("P", components, moleFractions, T=None, P=Plist[2])]

                T = [simResults[0][0], simResults[1][0], simResults[2][0]]
                P = [simResults[0][1], simResults[1][1], simResults[2][1]]

                betaGas = core.betaGas(T, P)

                T = [T[1]]
                P = [P[1]]
                simResults = [simResults[1]]
                freezingPoint = simResult[0][7]

                for i in range(len(T)):
                    if sum(inhibitorConcs) + sum(saltConcs) != 0:
                        TInhibited[i] = round(core.HuLeeSum(eqTemperature[0], saltConcs, inhibitorConcs, betaGas, freezingPoint), TDecimals)
                    else:
                        TInhibited[i] = eqTemperature[0]
                    if definedVariable == "T":
                        T[i] = round(simFunctions.tempConversion(tempUnit, T[i], False), TDecimals)
                    else:
                        T[i] = round(simFunctions.tempConversion(tempUnit, T[i][0], False), TDecimals)
                        eqTemperature[i] = simFunctions.tempConversion(tempUnit, eqTemperature[i], False)
                        eqTemperature[i] = round(eqTemperature[i], TDecimals)
                    if TInhibited[i] != None:
                        TInhibited[i] = round(simFunctions.tempConversion(tempUnit, TInhibited[i], False), TDecimals)
                    
                    eqPressure[i] = simFunctions.pressureConversion(pressureUnit, eqPressure[i], False)
                    eqPressure[i] = round(eqPressure[i], PDecimals) 

            endTime = time.time()
            st.text("Time to Complete Calculation: " + str(round(endTime-startTime, 3)) + " seconds")
            if len(T) > 1:
                st.text("per Data Point: " + str(round((endTime-startTime)/noPoints, 3)) + " seconds")
            calculated = True
        else:
            try:
                if sum(moleFractions[0]) != 1:
                    st.markdown(f":red[Sum of Mole Fractions is Not 1]")
            except:
                if sum(moleFractions) != 1:
                    st.markdown(f":red[Sum of Mole Fractions is Not 1]")
            if simFunctions.checkMaxConc(inhibitorConcs) != "":
                st.markdown(f":red[" + "Inhibitor(s) " + simFunctions.checkMaxConc(inhibitorConcs) + " Exceed(s) Maximum Concentration" + "]")
            if sum(inhibitorConcs)+sum(saltConcs) >= 100:
                st.markdown(f":red[Weight Percent of Inhibitors and Salts Exceeds 100%]")
    if calculated == True:
        if definedVariable == "T":
            P = eqPressure
            T = [round(val, TDecimals) for val in T]
        elif definedVariable == "P":
            T = eqTemperature
            P = [round(val, PDecimals) for val in P]

        if sum(inhibitorConcs) + sum(saltConcs) != 0:
            displayData = pd.DataFrame({'Temperature ('+tempUnit+')': T, 'Inhibited T ('+tempUnit+')': TInhibited, 'Pressure ('+pressureUnit+')': P, 'Eq. Structure': eqStructure}, [i for i in range(len(T))])
            data = pd.DataFrame(simFunctions.generateOutput(compounds, components, moleFractions, salts, saltConcs, inhibitors, inhibitorConcs, T, TInhibited, P, simResult, IDs, tempUnit, pressureUnit))
        else:
            displayData = pd.DataFrame({'Temperature ('+tempUnit+')': T, 'Pressure ('+pressureUnit+')': P, 'Eq. Structure': eqStructure}, [i for i in range(len(T))])
            data = pd.DataFrame(simFunctions.generateOutput(compounds, components, moleFractions, salts, saltConcs, inhibitors, inhibitorConcs, T, TInhibited, P, simResult, IDs, tempUnit, pressureUnit))
        st.dataframe(displayData, hide_index = True)
        c1, c2 = st.columns(2)
        with c1:
            st.download_button("Save Full Dataset", data=data.to_csv(index=False, header=False).encode('utf-8-sig'), file_name='data.csv', mime='text/csv')
        with c2:
            img = io.BytesIO()
            plt.savefig(img, format='png')
            downloadButton = st.download_button(label="Save Plot", data=img, file_name="Plot.png", mime="image/png")
        st.text("Full dataset includes temperature, inhibited temperature, pressure, structure, small and large cage occupancy, hydration number, hydrate density, and phase for each point")

        st.caption('NOTE: After selecting "Full Data Download" or "Download Plot", the page will appear to reset. If no changes are made to system parameters, just select "Calculate" again, and you can select other options as desired.')

    st.text("To reset this application, please refresh your page (F5 on Windows or Command+R on Mac)")

elif programType == "Minimum Concentration Calculation":
    tempUnit = st.radio("Temperature Unit", ["K", "°C", "°F"], horizontal=True)
    T = float(st.text_input("Fresh Water Equilibrium Temperature ("+tempUnit+"): ", value=round(simFunctions.tempConversion(tempUnit, 280, False),1)))
    T = [T+0.5, T, T-0.5]
    TDesired = float(st.text_input("Desired Inhibited Equilibrium Temperature ("+tempUnit+"): ", value=round(simFunctions.tempConversion(tempUnit, 275, False),1)))

    components = []
    moleFractions = []

    #Mole Fraction Input Table
    massFraction = st.toggle("Input Mass Fractions", False)
    c1, c2, c3 = st.columns([1, 0.9, 1.1], gap="small")
    with c1:
        if massFraction == False:
            compDf = pd.DataFrame([])
            compDf = pd.concat([compDf, pd.DataFrame([{'Component': componentList[0], 'Mole Fraction': 1.}])], ignore_index=True)
            for i in range(len(componentList)-1):
                compDf = pd.concat([compDf, pd.DataFrame([{'Component': componentList[i+1], 'Mole Fraction': 0.}])], ignore_index=True)

            inputCompDf = st.data_editor(compDf, hide_index=True, column_config={
                "Component": st.column_config.TextColumn("Component", disabled=True),
                "Mole Fraction": st.column_config.NumberColumn("Mole Fraction"),
            })

            moleFracInput = inputCompDf['Mole Fraction'].tolist()

            for i in range(len(componentList)):
                components.append(i + 1)
                moleFractions.append(round(moleFracInput[i],4))
        else:
            massFractions = []
            compDf = pd.DataFrame([])
            compDf = pd.concat([compDf, pd.DataFrame([{'Component': componentList[0], 'Mass Fraction': 1.}])], ignore_index=True)
            for i in range(len(componentList)-1):
                compDf = pd.concat([compDf, pd.DataFrame([{'Component': componentList[i+1], 'Mass Fraction': 0.}])], ignore_index=True)

            inputCompDf = st.data_editor(compDf, hide_index=True, column_config={
                "Component": st.column_config.TextColumn("Component", disabled=True),
                "Mass Fraction": st.column_config.NumberColumn("Mass Fraction"),
            })

            massFracInput = inputCompDf['Mass Fraction'].tolist()

            for i in range(len(componentList)):
                components.append(i + 1)
                massFractions.append(round(massFracInput[i],4))

    normalizeFracs = st.toggle("Normalize Mole Fractions", False)

    with c2:
        if massFraction == True:
            moleFractions = simFunctions.massToMolFrac(compounds, massFractions)

            normDf = pd.DataFrame([])
            for i in range(len(componentList)):
                normDf = pd.concat([normDf, pd.DataFrame([{'Component': componentList[i], 'Mole Fraction': moleFractions[i]}])], ignore_index=True)

            inputCompDf = st.dataframe(normDf, hide_index=True)

    with c3:
        if normalizeFracs == True:
            nonNormalSum = sum(moleFractions)

            for i in range(len(moleFractions)):
                moleFractions[i]/=nonNormalSum

            normDf = pd.DataFrame([])
            for i in range(len(componentList)):
                normDf = pd.concat([normDf, pd.DataFrame([{'Component': componentList[i], 'Normalized Mole Frac.': moleFractions[i]}])], ignore_index=True)

            inputCompDf = st.dataframe(normDf, hide_index=True)

    nonZeroFracs = numpy.nonzero(moleFractions)[0]
    moleFractions = [element for index, element in enumerate(moleFractions) if index in nonZeroFracs]
    components = [element for index, element in enumerate(components) if index in nonZeroFracs]

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
        simResults = [equilibrium("T", components, moleFractions, T=T[0], P=None),
             equilibrium("T", components, moleFractions, T=T[1], P=None),
             equilibrium("T", components, moleFractions, T=T[2], P=None)]
        
        P = [simResults[0][1], simResults[1][1], simResults[2][1]]

        freezingPoints = [simResults[0][7], simResults[1][7], simResults[2][7]]
        freezingPoint = min(freezingPoints, key=lambda x: abs(x - 273.15))

        betaGas = core.betaGas(T, P)

        TDesired = simFunctions.tempConversion(tempUnit, TDesired, False)
        if inhibitor != "salt":
            conc = simFunctions.getConcentration(T[1], TDesired, inhibitor, salt, betaGas, freezingPoint)
            st.text("Minimum Concentration of " + str(inhibitorList[inhibitor]) + ": " + str(round(conc,1)) + "% w/w")
        else:
            conc = simFunctions.getConcentration(T[1], TDesired, inhibitor, salt, betaGas, freezingPoint)
            st.text("Minimum Concentration of " + str(saltList[salt]) + ": " + str(round(conc,1)) + "% w/w")

st.caption("Disclaimer: The model and predictions have been tested and verified to be accurate based on extensive comparison with available literature data. However, this web app is provided ""as is"" and ""as available"" without any warranties of any kind, either express or implied, including, but not limited to, implied warranties of use, merchantability, fitness for a particular purpose, and non-infringement.")

st.markdown('''Questions, suggestions, or bug reports? Please send an email to asum@mines.edu''')

st.markdown('''To cite this software, please cite the following article: [[doi: 10.1016/j.softx.2025.102422]](https://www.sciencedirect.com/science/article/pii/S2352711025003887)''')

st.header('Credits')
st.markdown('''
            Developed by Karsten Kunneman in collaboration with Prof. Amadeu K. Sum at the Colorado School of Mines
            \nHydrate model: Klauda-Sandler fugacity model [[doi: 10.1021/ie000322b]](https://doi.org/10.1021/ie000322b)
            \nInhibition model: HLS correlation [[doi: 10.1002/aic.16369]](https://doi.org/10.1002/aic.16369)
            \nThis site is created with Streamlit''')
st.markdown(f'''<a href="https://github.com/karstenkunneman/P2F-HydrateCalc">Github Repo</a>''', unsafe_allow_html=True)
