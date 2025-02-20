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
    componentList.append(str(IDs[i]) + ". " + compounds[i])

st.title('Gas Hydrate Equilibrium Calculator')
st.caption('Version 2025-02-27')

programType = st.radio("Calculation Type", ["Equilibrium Calculator", "Minimum Concentration Calculator"], horizontal=True)

if programType == "Equilibrium Calculator":
    noComponents = st.number_input('Number of Components', 1, None, 1, 1)
    components = []
    moleFractions = []
    c1, c2 = st.columns(2)
    with c1:
        for i in range(noComponents):
            components += [int(st.selectbox('Component ' + str(i+1), componentList)[0])]
    with c2:
        if noComponents > 1:
            for i in range(noComponents):
                moleFractions += [float(st.text_input('Component ' + str(i+1) + " Mole Fraction", 1.00))]
        else:
            moleFractions = [1]

    csvGuesses = st.file_uploader("Upload Temperatures and Guess Pressures", ['csv'])
    csvTemplate = st.download_button("Guess File Template", open("Input Template.csv", encoding='utf-8'), file_name="Input Template.csv")

    if csvGuesses == None:
        calculateRange = st.toggle('Calculate Range of Temperatures', False)
        if calculateRange == True:
            minTemp = float(st.text_input('Minimum Temperature (K)', 190.15))
            maxTemp = float(st.text_input('Maximum Temperature (K)', 302.1))
            noPoints = st.number_input('Number of Points', 1, None, 4, 1)
            T = numpy.arange(maxTemp, minTemp-(maxTemp-minTemp)/noPoints, -1*(maxTemp-minTemp)/(noPoints-1))

            minGuessPressure = float(st.text_input("Minimum Guess Pressure (MPa): ", 0.089))*1E6 #Pressure in Pa
            maxGuessPressure = float(st.text_input("Maximum Guess Pressure (MPa): ", 74.291))*1E6 #Pressure in Pa
            logP = numpy.arange(math.log(maxGuessPressure), math.log(minGuessPressure)-(math.log(maxGuessPressure)-math.log(minGuessPressure))/noPoints, -1*(math.log(maxGuessPressure)-math.log(minGuessPressure))/(noPoints-1))
            P = [0 for i in range(len(T))]
            for i in range(len(logP)):
                T[i] = round(T[i], 2)
                P[i] = math.exp(logP[i])
        else:
            T = [float(st.text_input('Temperature (K)', 278.1))]
            P = [float(st.text_input('Guess Pressure (MPA)', 4.249))*1E6]
    else:
        guessFile = numpy.genfromtxt(csvGuesses, delimiter=',', skip_header=1)
        T = guessFile[:,0]
        P = guessFile[:,1]*1E6
        noPoints = len(T)
        if noPoints > 1:
            calculateRange = True
        else:
            calculateRange = False

    salts, inhibitors = simFunctions.getInhibitors()
    saltConcs = [0 for i in range(len(salts))]
    inhibitorConcs = [0 for i in range(len(inhibitors))]
    freshWater = st.toggle('Fresh Water', True)
    if freshWater == False:
        betaGas = float(st.text_input("betaGas: ", 9.120E-04))
        c1, c2 = st.columns(2)
        with c1:
            for i in range(len(salts)):
                saltConcs[i] = float(st.text_input(str(salts[i] + " Concentration (wt%): "), 0))
        with c2:
            for i in range(len(inhibitors)):
                inhibitorConcs[i] = float(st.text_input(str(inhibitors[i] + " Concentration (wt%): "), 0))

    if st.button("Calculate"):
        startTime = time.time()
        if calculateRange == True:
            eqPressure = [0 for i in range(len(T))]
            eqStructure = [0 for i in range(len(T))]
            eqFractions = [0 for i in range(len(T))]
            TInhibited = [0 for i in range(len(T))]
            progressBar = st.progress(0, str(0) + "/" + str(len(T)))
            for i in range(len(T)):
                simResult = simFunctions.equilibriumPressure(T[i], P[i], components, moleFractions, saltConcs, inhibitorConcs)
                eqPressure[i] = simResult[0]/1E6 #In MPa
                eqStructure[i] = simResult[1]
                eqFractions[i] = simResult[2]
                if freshWater == False:
                    TInhibited[i] = simFunctions.HuLeeSum(T[i], saltConcs, inhibitorConcs, betaGas)
                else:
                    TInhibited[i] = T[i]
                with progressBar:
                    st.progress((i+1)/len(T), str(i+1) + "/" + str(len(T)))
            data = pd.DataFrame({'Temp (K)': T, 'Inhibited Temp (K)': TInhibited, 'Eq. Pressure (MPa)': eqPressure, 'Eq. Structure': eqStructure, 'Occupancy Fractions (Small Cage, Large Cage)': eqFractions}, [i for i in range(len(T))])
            fig = plt.figure()
            plt.plot(T, eqPressure, '-', label='Fresh Water')
            if freshWater == False:
                plt.plot(TInhibited, eqPressure, '--', label='Inhibited System')
                plt.legend()
            plt.yscale("log")
            plt.xlabel("Temperature (K)")
            plt.ylabel("Pressure (MPa)")
            st.pyplot(fig)
        else:
            simResult = simFunctions.equilibriumPressure(T[0], P[0], components, moleFractions, saltConcs, inhibitorConcs)
            eqPressure = simResult[0]/1E6 #In MPa
            eqStructure = simResult[1]
            eqFractions = simResult[2]
            TInhibited = simFunctions.HuLeeSum(T[0], saltConcs, inhibitorConcs, betaGas)
            data = pd.DataFrame({'Temp (K)': T[0], 'Inhibited Temp (K)': TInhibited, 'Eq. Pressure (MPa)': eqPressure, 'Eq. Structure': eqStructure, 'Occupancy Fractions (Small Cage, Large Cage)': [eqFractions]})
        
        st.dataframe(data, hide_index = True)
        endTime = time.time()
        st.text("Time to Complete Calculation: " + str(round(endTime-startTime, 3)) + " seconds")
        if len(T) > 1:
            st.text("per Data Point: " + str(round((endTime-startTime)/noPoints, 3)) + " seconds")
else:
    T = float(st.text_input("Fresh Water Equilibrium Temperature (K): ", value="280"))
    TDesired = float(st.text_input("Desired Inhibited Equilibrium Temperature (K): ", value="275"))
    betaGas = float(st.text_input("betaGas: ", 9.120E-04))
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