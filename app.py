import streamlit as st
import pandas as pd
import numpy as np

componentList = ['1. Methane', '2. Ethane', '3. Propane', '4. Isobutane', '5. Cyclopropane', '6. Hydrogen Sulfide', '7. Nitrogen', 'Carbon Dioxide']

st.title('Gas Hydrate Equilibrium Calculator')
st.caption('Version 2024-4-1')

temperature = float(st.text_input('Temperature (K)', 298))
noComponents = st.number_input('Number of Components', 1, None, 1, 1)
components = []
moleFractions = []
for i in range(noComponents):
    components += [int(st.selectbox('Component ' + str(i+1), componentList)[0])]
    moleFractions += [float(st.text_input('Component ' + str(i+1) + " Mole Fraction", 1.00))]



st.header('Credits')