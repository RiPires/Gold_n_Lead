#RiP
import csv
from matplotlib.pylab import *
import matplotlib.pyplot as plt
import sys

def VdG_dat2Lists(File):
    """
    Converts .dat files from VdG RBS into yield and channel lists
    INPUTS:
        "FILENAME.dat"
    OUTPUTS:
        Yield and Channel lists
    HOW TO USE:
        MyYield, MyChannel = VdG_dat2Lists("MyFile.mca")
    """
    with open(File, 'r') as file:
        reader = csv.reader(file, delimiter="\n", skipinitialspace=True)
        data = list(reader)
    ch = []
    y = []
    aux = []
    for i in range(128):
        aux.append(data[i][0].split())
    #print(aux)
    for i in range(len(aux)):
        for k in range(8):
            ch.append(int(8*(i)+k+1)) ## axes in channel
            y.append(float(aux[i][k]))
    
    return y, ch

#Yield, Channel = MCA2Lists('RBS_29Agosto/0822/RBS1run2.dat')


def VdG_dat2SIMNRAfile(File):
    """
    Converts .dat files from VdG RBS into yield single list file to be read with SIMNRA
    INPUTS:
        "FILENAME.dat"
    OUTPUTS:
        "FILENAME_SIMNRA.dat
    HOW TO USE:
        VdG_dat2SIMNRAfile(FILENAME)
    """
    Output = open('output_SIMNRA.dat', 'w')
    with open(File, 'r') as file:
        reader = csv.reader(file, delimiter="\n", skipinitialspace=True)
        data = list(reader)
    y = []
    aux = []
    for i in range(128):
        aux.append(data[i][0].split())
    #print(aux)
    for i in range(len(aux)):
        for k in range(8):
            Output.write(str(int(aux[i][k]))+'\n')
    Output.close()
    
    return 'Done'

VdG_dat2SIMNRAfile('0922/RBS1run14.dat')