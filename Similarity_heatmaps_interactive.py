# -*- coding: utf-8 -*-
"""
Created on Tue Jun 15 19:23:28 2021

@author: Lucas
"""

import streamlit as st
import pandas as pd
import base64
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs # para calcular Tanimoto similarity
import plotly.graph_objects as go
from pathlib import Path
import seaborn as sns; sns.set_theme()
import matplotlib.pyplot as plt


#---------------------------------#
# Page layout
## Page expands to full width
st.set_page_config(page_title='LIDEB Molecular similarity heatmap',
    layout='wide')

######
# Funcion para poner una imagen    
def img_to_bytes(img_path):
    img_bytes = Path(img_path).read_bytes()
    encoded = base64.b64encode(img_bytes).decode()
    return encoded

header_html = "<img src='data:image/png;base64,{}' class='img-fluid'>".format(img_to_bytes("cropped-header.png"))
st.markdown(header_html, unsafe_allow_html=True,)
######


#---------------------------------#
st.write("""
# LIDeB Tools

**WebApp to generate heatmaps based in the molecular similarity**

*Morgan Fingerprints and Tanimoto Coefficient* are used for the generation of **Similarity Heatmaps**.
""")


#---------------------------------#
# Sidebar - Collects user input features into dataframe
st.sidebar.header('Upload your datasets (Smiles)')

uploaded_file_1 = st.sidebar.file_uploader("Upload your Dataset 1 in a TXT file", type=["txt"])
uploaded_file_2 = st.sidebar.file_uploader("Upload your Dataset 2 in a TXT file", type=["txt"])

# Sidebar - Specify parameter settings
st.sidebar.header('ECFP Radio')
split_size = st.sidebar.slider('Morgan fingerprint Radio', 2, 4, 2, 1)

st.sidebar.subheader('ECFP LENGH')
parameter_n_estimators = st.sidebar.slider('Set the fingerprint lenght', 512, 2048, 1024, 512)

st.sidebar.header('Type of Plot')
type_plot = st.sidebar.checkbox('Interactive Plot')
if type_plot == True:
    plotly_color = st.sidebar.selectbox("Select the heatmap color", 
                         ('Blackbody','Bluered','Blues','Earth','Electric','Greens',
                          'Greys','Hot','Jet','Picnic','Portland','Rainbow','RdBu','Reds','Viridis','YlGnBu','YlOrRd'),
                         15)
else:
    sns_color = st.sidebar.selectbox("Select the heatmap color", ("rocket", "mako", "flare","crest","magma","viridis"),5)



#---------------------------------#
# Main panel

# Displays the dataset
st.subheader('Dataset')

#%%

# GENERAMOS LAS FINGERPRINTS DE CADA DATASET POR SEPARADO

def similarity(df_1,df_2):
    df_1 = df_1[0].tolist()
    df_2 = df_2[0].tolist()
    # df_2 = df_1.copy()
    lenght_dataset_1 = len(df_1)
    lenght_dataset_2 = len(df_2)
    st.markdown('Dataset 1 have: ' + str(lenght_dataset_1) + " molecules")
    st.markdown('Dataset 2 have: ' + str(lenght_dataset_2) + " molecules")
    
    fps_1 = []    
    for m in df_1:
        mol = Chem.MolFromSmiles(m)
        fp_1 = AllChem.GetMorganFingerprintAsBitVect(mol,split_size,nBits=parameter_n_estimators,useFeatures=True)
        fps_1.append(fp_1)
    
    fps_2 = []
    for m1 in df_2:
        mole = Chem.MolFromSmiles(m1)
        fp_2 = AllChem.GetMorganFingerprintAsBitVect(mole,split_size,nBits=parameter_n_estimators,useFeatures=True)
        fps_2.append(fp_2)
    
    # COMPARAMOS LAS FINGERPRINTS POR TANIMOTO Y GENERAMOS LA MATRIZ
    
    matriz_tanimoto = []
    for finger1 in fps_1:
        lista = []
        for finger2 in fps_2:
            tan_sim_ac_in=DataStructs.TanimotoSimilarity(finger1, finger2)
            lista.append(tan_sim_ac_in)
        matriz_tanimoto.append(lista)
    
    df_ok = pd.DataFrame(matriz_tanimoto)
    filas = list(range(1,len(fps_1)+1,1))
    columnas = list(range(1,len(fps_2)+1,1))
    df_ok.index = filas
    df_ok.columns = columnas
    return df_ok


def heatmap(df_ok):
    #-----Plot-----#
    if type_plot == True:
        # color = "YlGnBu"
        fig = go.Figure(go.Heatmap(z=df_ok,x0=1,dx=1, y0=1,dy=1, hoverongaps = False,showscale=True, colorscale=plotly_color,zmax=1,zmin=0))
        st.plotly_chart(fig)
    else:
        df_inv = df_ok.T
        ax = sns.heatmap(df_inv, xticklabels=False, yticklabels=False,linewidths=.3,cmap=sns_color)
        plt.xlabel("Dataset 2")
        plt.ylabel("Dataset 1")
        st.set_option('deprecation.showPyplotGlobalUse', False)
        st.pyplot()
        return ax

# PARA GENERAR LA MATRIX Y EL HEATMAP

# ---------------------------------#

if uploaded_file_1 is not None and uploaded_file_2 is not None:
    df_1 = pd.read_csv(uploaded_file_1,sep="\t",header=None)
    df_2 = pd.read_csv(uploaded_file_2,sep="\t",header=None)
    df_ok = similarity(df_1,df_2)
    plot = heatmap(df_ok)
    st.markdown("You can download the heatmap by Right Click in the image and then **'save image as'** :blush: ")
    
# Example file
else:
    st.info('Awaiting for TXT file to be uploaded.')
    if st.button('Press to use Example Dataset'):
        df_1 = pd.read_csv("molecules_1.txt",sep="\t",header=None)
        df_2 = pd.read_csv("molecules_1.txt",sep="\t",header=None)
        df_ok = similarity(df_1,df_2)    
        plot = heatmap(df_ok)
        st.markdown('A dataset of **40 smiles** has been used as the example.')
        st.markdown("You can download the heatmap by Right Click in the image and then **'save image as'** :blush: ")

