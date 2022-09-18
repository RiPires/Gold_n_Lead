#RiP
from matplotlib.pylab import *
import matplotlib.pyplot as plt
import csv


#########################################
###   Plot data funntion definition   ###
#########################################
def PlotData(File):
    """
    Plots yield vs channel data from our .mca file

    INPUTS: "FileName.mca"
    OUTPUTS: yield vs channel plot
    """

    with open(File, 'r') as file:
        reader = csv.reader(file, delimiter="\n", skipinitialspace=True)
        data = list(reader)
    x = []
    y = []
    aux = []
    for i in range(12,len(data)-1):
        aux.append(data[i][0].split())
    for i in range(len(aux)):
        x.append(float(0.004482*i-0.105144)) ## energy calibraion corrected
        #x.append(float(i)) ## axes in channel
        y.append(float(aux[i][0]))

    fig, ax = plt.subplots()
    ax.plot(x,y,'.', color ='xkcd:black', label=(str(File)))
    #ax.semilogy(x,y,'^', color ='xkcd:purple', label=(str(File)))
    legend = ax.legend(loc="upper right",ncol=2, shadow=False,fancybox=True,framealpha = 0.0,fontsize=20)
    legend.get_frame().set_facecolor('#DAEBF2')
    tick_params(axis='both', which='major', labelsize=22)
    xlabel('Energy (MeV)',fontsize=22)
    ylabel('Counts', fontsize=22)
    grid()
    show()
#########################################
#########################################


#############################################
###   Plot RBS .dat funntion definition   ###
#############################################
def VdG_Plot(File):
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
    
    fig, ax = plt.subplots()
    ax.plot(ch,y,'.', color ='xkcd:black', label=(str(File)))
    #ax.semilogy(ch,y,'*', color ='xkcd:purple', label=(str(File)))
    legend = ax.legend(loc="upper right",ncol=2, shadow=False,fancybox=True,framealpha = 0.0,fontsize=20)
    legend.get_frame().set_facecolor('#DAEBF2')
    tick_params(axis='both', which='major', labelsize=22)
    xlabel('Energy (keV)',fontsize=22)
    ylabel('Yield', fontsize=22)
    grid()
    show()
#############################################
#############################################

###                 Plots                     ###
VdG_Plot('20220907/ERDrun01.dat') # V-Ta-Nb    
VdG_Plot('20220907/ERDrun03.dat') # Yb-Ge-Si-O  
#################################################