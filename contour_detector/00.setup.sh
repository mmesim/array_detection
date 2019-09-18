#!/bin/bash
# Create anaconda environment and install packages ############ 
# for frequency detector algorithm                 ############
################################   M.M. 04/2019    ############ 
#Create environment
yes |  conda create --name freqtor python=2.7.11
#Qddd conda-forge channel
conda config --append channels conda-forge
#Load environment
source  $HOME/anaconda3/bin/activate freqtor 
#Obspy
yes | conda install obspy=1.0.2
#Pandas
yes | conda install pandas=0.19.2
#Geopy
yes | conda install geopy=1.11.0 
#Basemap
yes | conda install basemap=1.0.7
#Ipython
yes | conda install Ipython=4.1.2
#Ipykernel
yes | conda install Ipykernel=4.3.1
#Matplotlib
yes | conda install matplotlib=1.5.2
