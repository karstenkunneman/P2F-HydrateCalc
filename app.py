import simFunctions
import streamlit as st
import pandas as pd
import numpy
import matplotlib.pyplot as plt
import time

IDs, compounds = simFunctions.getComponents()
componentList = []
for i in range(len(IDs)):
    componentList.append(str(IDs[i]) + ". " + compounds[i])

st.title('Gas Hydrate Equilibrium Calculator')
st.caption('Version 2025-02-13')

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
        minTemp = float(st.text_input('Minimum Temperature (K)', 148.8))
        maxTemp = float(st.text_input('Maximum Temperature (K)', 270.9))
        noPoints = st.number_input('Number of Points', 1, None, 8, 1)
        T = numpy.arange(minTemp, maxTemp+(maxTemp-minTemp)/noPoints, (maxTemp-minTemp)/(noPoints-1))

        minGuessPressure = float(st.text_input("Minimum Guess Pressure (MPa): ", 0.056))*1E6 #Pressure in Pa
        maxGuessPressure = float(st.text_input("Maximum Guess Pressure (MPa): ", 2.39))*1E6 #Pressure in Pa
        P = numpy.arange(maxGuessPressure, minGuessPressure-(maxGuessPressure-minGuessPressure)/noPoints, -1*(maxGuessPressure-minGuessPressure)/(noPoints-1))
    else:
        T = [float(st.text_input('Temperature (K)', 298))]
        P = [float(st.text_input('Guess Pressure (MPA)', 1.00))*1E6]
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
        progressBar = st.progress(0, str(0) + "/" + str(len(T)))
        for i in range(len(T)):
            simResult = simFunctions.equilibriumPressure(T[i], P[i], components, moleFractions, saltConcs, inhibitorConcs)
            eqPressure[i] = simResult[0]/1E6 #In MPa
            eqStructure[i] = simResult[1]
            eqFractions[i] = simResult[2]
            with progressBar:
                st.progress((i+1)/len(T), str(i+1) + "/" + str(len(T)))
        data = pd.DataFrame({'Temp (K)': T, 'Eq. Pressure (MPa)': eqPressure, 'Eq. Structure': eqStructure, 'Occupancy Fractions': eqFractions}, [i for i in range(len(T))])
        fig = plt.figure()
        plt.plot(T, eqPressure, '-ok')
        plt.yscale("log")
        plt.title("Equilibrium Predictions")
        plt.xlabel("Temperature (K)")
        plt.ylabel("Pressure (MPa)")
        st.pyplot(fig)
    else:
        simResult = simFunctions.equilibriumPressure(T[i], P[i], components, moleFractions, saltConcs, inhibitorConcs)
        eqPressure = simResult[0]/1E6 #In MPa
        eqStructure = simResult[1]
        eqFractions = simResult[2]
        data = pd.DataFrame({'Temp (K)': T[0], 'Eq. Pressure (MPa)': eqPressure, 'Eq. Structure': eqStructure, 'Occupancy Fractions': [eqFractions]})
    
    st.dataframe(data, hide_index = True)
    endTime = time.time()
    st.text("Time to Complete Calculation: " + str(round(endTime-startTime, 3)) + " seconds")
    if len(T) > 1:
        st.text("per Data Point: " + str(round((endTime-startTime)/noPoints, 3)) + " seconds")    

st.header('Credits')
st.markdown('''
            Created by Karsten Kunneman and Dr. Amadeu K Sum at the Colorado School of Mines  
            Based on "A Fugacity Model for Gas Hydrate Phase Equilibria" by Jeffery B. Klauda and Stanley I. Sandler
            This site created with Streamlit''')
st.markdown(f'''<a href="https://github.com/karstenkunneman/Gas-Hydrate-Equilibrium-Calculator">Github Repo</a>''', unsafe_allow_html=True)