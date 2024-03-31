import simulation
import streamlit as st
import pandas as pd
import numpy
import matplotlib.pyplot as plt

componentList = ['1. Methane', '2. Ethane', '3. Propane', '4. Isobutane', '5. Cyclopropane', '6. Hydrogen Sulfide', '7. Nitrogen', 'Carbon Dioxide']

st.title('Gas Hydrate Equilibrium Calculator')
st.caption('Version 2024-4-1')

noComponents = st.number_input('Number of Components', 1, None, 1, 1)
components = []
moleFractions = []
for i in range(noComponents):
    components += [int(st.selectbox('Component ' + str(i+1), componentList)[0])]
    moleFractions += [float(st.text_input('Component ' + str(i+1) + " Mole Fraction", 1.00))]

guessPressure = float(st.text_input('Guess Pressure (MPA)', 1.00))*1E6

calculateRange = st.toggle('Calculate Range of Temperatures', False)

if calculateRange == True:
    minTemp = float(st.text_input('Minimum Temperature (K)', 275))
    maxTemp = float(st.text_input('Maximum Temperature (K)', 300))
    noPoints = st.number_input('Number of Points', 1, None, 8, 1)
    T = numpy.arange(minTemp, maxTemp+(maxTemp-minTemp)/noPoints, (maxTemp-minTemp)/(noPoints-1))
else:
    T = [float(st.text_input('Temperature (K)', 298))]
    
if st.button("Calculate"):
    if calculateRange == True:
        eqPressure = [0 for i in range(len(T))]
        for i in range(len(T)):
            eqPressure[i] = simulation.equilibriumPressure(T[i], guessPressure, components, moleFractions)/1E6 #In MPa
        fig = plt.figure()
        plt.plot(T, eqPressure, '-ok')
        plt.yscale("log")
        plt.title("Equilibrium Predictions")
        plt.xlabel("Temperature (K)")
        plt.ylabel("Pressure (MPa)")
        st.pyplot(fig)
    else:
        eqPressure = [simulation.equilibriumPressure(T[0], guessPressure, components, moleFractions)/1E6]
    data = pd.DataFrame({'Temp (K)': T, 'Eq. Pressure (MPa)': eqPressure}, [i for i in range(len(T))])
    st.dataframe(data, hide_index = True)

st.header('Credits')
st.markdown('''
            Created by Karsten Kunneman and Dr. Amadeu K Sum at the Colorado School of Mines  
            Based on "A Fugacity Model for Gas Hydrate Phase Equilibria" by Jeffery B. Klauda and Stanley I. Sandler  
            Compound data obrained from NIST  
            Equations used from "PRSV: An Improved Peng- Robinson Equation of State for Pure Compounds and Mixtures" by R. Stryjek and J. H. Vera  
            Binary interaction parameters obtained from "Vapor-liquid equilibria for mixtures of low boiling substances. Pt. 2. Ternary systems" by H. Knapp, S. Zeck, and R. Langhorst  
            This site created with Streamlit
            ''')