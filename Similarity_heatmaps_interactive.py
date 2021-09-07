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

from PIL import Image
image = Image.open('cropped-header-heatmaps.png')
st.image(image)

st.markdown("![Twitter Follow](https://img.shields.io/twitter/follow/LIDeB_UNLP?style=social)")
st.subheader(":pushpin:" "About Us")
st.markdown("We are a team interested to develop new tools cheminformatics for using in areas of computer-assisted drug design and machine learning in drug discovery. We belong to Laboratory of Bioactive Research and Development (LIDeB) of the National University of La Plata, Argentina. Our research group is focused on rational designing of new drug candidates and drug repurposing for epilepsy and neglected tropical diseases such as Chagas disease, leishmaniasis, malaria. Another important goal of our group is the development and caracterization of nanocarriers. The work developed by our group has resulted in publications in international indexed journals, abstracts, congress and awards in national and international of Medicinal and Computational Chemistry.")
st.markdown(":computer:""**Web Site** " "<https://lideb.biol.unlp.edu.ar>")


#---------------------------------#
st.write("""
# LIDeB Tools - Similarity Heatmaps

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

st.sidebar.title(":speech_balloon: Contact Us")
st.sidebar.info(
"""
If you are looking to contact us, please
[:e-mail:](mailto:lideb@biol.unlp.edu.ar) or [T](https://twitter.com/LIDeB_UNLP)
""")

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
        fig.update_xaxes(title_text='Dataset 2')
        fig.update_yaxes(title_text='Dataset 1')
        fig.update_layout(margin = dict(t=60,r=20,b=20,l=20),
        width = 800, height = 800,
        autosize = False )
        st.plotly_chart(fig)
    else:
        df_ok.sort_index(axis=0, ascending=False,inplace=True)
        ax = sns.heatmap(df_ok, xticklabels=False, yticklabels=False,cmap=sns_color)
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

        
#Footer edit

footer="""<style>
a:link , a:visited{
color: blue;
background-color: transparent;
text-decoration: underline;
}
a:hover,  a:active {
color: red;
background-color: transparent;
text-decoration: underline;
}
.footer {
position: fixed;
left: 0;
bottom: 0;
width: 100%;
background-color: white;
color: black;
text-align: center;
}
</style>
<div class="footer">
<p>Made in  🐍 and <img style='display: ; ' href="https://streamlit.io" src="https://i.imgur.com/iIOA6kU.png" target="_blank"></img> Developed with ❤️ by <a style='display: ; text-align: center' href="https://twitter.com/capigol" target="_blank">Lucas Alberca</a> for <a style='display:; text-align: center;' href="https://lideb.biol.unlp.edu.ar/" target="_blank">LIDeB</a></p>
</div>
"""
st.markdown(footer,unsafe_allow_html=True)


