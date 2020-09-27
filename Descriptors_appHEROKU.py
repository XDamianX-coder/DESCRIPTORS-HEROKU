######################
# Import libraries
######################
import numpy as np
import pandas as pd
import streamlit as st
import time
from PIL import Image
import matplotlib.pyplot as plt
import pubchempy as pcp

######################
# Custom function
######################
## Calculate molecular descriptors
from typing import Any



def cid_to_smiles(smiles):
    baseData_1=[]
    p= pcp.get_compounds(smiles, 'smiles')
    baseData_1.append(p)
    baseData_1 = np.arange(1, 1)
    i = 0
    row = np.array([p])
    if (i == 0):
        baseData_1 = row
    else:
        baseData_1 = np.vstack([baseData_1, row])
    i = i + 1
    columnNames = ["Compound CID"]
    cid = pd.DataFrame(data=baseData_1, columns=columnNames)
    return cid

def molecular_formula(number):
    baseData_2 = []
    g = pcp.Compound.from_cid(number)
    #print(g.molecular_formula)
    baseData_2.append(g.molecular_formula)
    baseData_2 = np.arange(1, 1)
    i = 0
    row = np.array([g.molecular_formula])
    if (i == 0):
        baseData_2 = row
    else:
        baseData_2 = np.vstack([baseData_2, row])
    i = i + 1
    columnNames = ["Molecular Formula"]
    MF = pd.DataFrame(data=baseData_2, columns=columnNames)
    return MF

def iupac_name(number):
    baseData_3 = []
    u = pcp.Compound.from_cid(number)
    baseData_3.append(u.iupac_name)
    baseData_3 = np.arange(1, 1)
    i = 0
    row = np.array([u.iupac_name])
    if (i == 0):
        baseData_3 = row
    else:
        baseData_3 = np.vstack([baseData_3, row])
    i = i + 1
    columnNames = ["IUPAC name"]
    IUPAC = pd.DataFrame(data=baseData_3, columns=columnNames)

    return IUPAC


def search_by_cpnd_name(name):

    result = pcp.get_compounds(name,'name')
    baseData_4 = []
    baseData_4.append(result)
    baseData_4 = np.arange(1, 1)
    i = 0
    row = np.array([result])
    if (i == 0):
        baseData_4 = row
    else:
        baseData_4 = np.vstack([baseData_4, row])
    i = i + 1
    columnNames = ["List of compounds"]
    List_cpnds = pd.DataFrame(data=baseData_4, columns=columnNames)

    return List_cpnds



def compound_properties(CIDs):

    jk = pcp.get_compounds(CIDs, namespace=u'cid', searchtype=None, as_dataframe=True)
    return jk


######################
# Page Title
######################


st.write("""
# Molecular Descriptors Web App
This app calculate the **Molecular Descriptors** of compounds!\n
Author: **Damian Nowak**
Github repo: \n
https://github.com/XDamianX-coder/Molecular_descriptors_app \n
Molecular_descriptors_app contains a bunch of new features \n
which can't be present in online version. \n
Contact: dam.now22@protonmail.com \n
You can also check my other applications.\n
The random number generator in given range.\n
https://github.com/XDamianX-coder/Random_numbers \n
Calculation of the surface areas of figures and volume of solids. \n
https://github.com/XDamianX-coder/Math_areas_and_volumes
""")

my_bar = st.progress(0)
for percent_complete in range(100):
    time.sleep(0.03)
    my_bar.progress(percent_complete + 1)


######################
# Input molecules (Side Panel)
######################

st.sidebar.header('User Input Features')

## Read SMILES input



SMILES_input_1 = "O=C1CN=C(C2=CC=CC=C2)C2=C(N1)C=CC=C2"
SMILES_1 = st.sidebar.text_area("SMILES input to get CID", SMILES_input_1)
SMILES_1 = SMILES_1.split('\n')

st.header('Input SMILES to get the compound cid.')
SMILES_1
A = cid_to_smiles(SMILES_1)
A





CID_input = "45266826"
CID = st.sidebar.text_area("CID input to get molecular formula.", CID_input)
CID = CID.split('\n')

st.header('Input CID to get molecular formula.')
CID
B = molecular_formula(CID)
B

CID_input_1 = "3007"
CID_1 = st.sidebar.text_area("CID input to get IUPAC name.", CID_input_1)
CID_1 = CID_1.split('\n')

st.header('Input CID to get IUPAC name of compound.')
CID_1
C = iupac_name(CID_1)
C

Name_input = "Metamphetamine"
Name_1 = st.sidebar.text_area("Name input to CID number.", Name_input)
Name_1 = Name_1.split('\n')

st.header('Input name to CID number.')
Name_1
D = search_by_cpnd_name(Name_1)
D


CID_input_1 ="3007, 134664, 76175, 2476, 7843, 107786, 2826, 1615, 702, 10413, 3592, 3821, 3981, 45266826, 4095, 1206, 2978"
CIDs = st.sidebar.text_area("CIDs input to get compound properties.", CID_input_1)
CIDs = CIDs.split(', ')
st.header('Input CIDS to get compound properties.')
CIDs
E = compound_properties(CIDs)
E

time.sleep(0.1)
plot_1 = compound_properties(CIDs)
fig, ax = plt.subplots()
ax.hist(plot_1['exact_mass'], bins=100)
st.header('Distribution of exact mass.')
st.pyplot(fig)

time.sleep(0.1)
plot_1 = compound_properties(CIDs)
fig, ax = plt.subplots()
ax.hist(plot_1['tpsa'], bins=100)
st.header('Distribution of tpsa.')
st.pyplot(fig)


time.sleep(0.1)
plot_1 = compound_properties(CIDs)
fig, ax = plt.subplots()
ax.hist(plot_1['complexity'], bins=100)
st.header('Distribution of complexity.')
st.pyplot(fig)

time.sleep(0.1)
plot_1 = compound_properties(CIDs)
fig, ax = plt.subplots()
ax.hist(plot_1['heavy_atom_count'], bins=100)
st.header('Distribution of number of heavy atoms.')
st.pyplot(fig)
