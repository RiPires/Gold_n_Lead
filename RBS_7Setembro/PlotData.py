#RiP
from matplotlib.pylab import *
import matplotlib.pyplot as plt
import csv

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

#######################################################
###   Plot two RBS .dat files funntion definition   ###
#######################################################
def VdG_PlotBoth(File1, File2):
    """
    Converts two .dat files from VdG RBS into yield and channel lists
    and plots both
    INPUTS:
        "FILENAME1.dat",'FILENAME2.dat'
    OUTPUTS:
        Plots
    """
    ## Creates channel yield lists for File1 and
    # converts channel to energy
    with open(File1, 'r') as file:
        reader = csv.reader(file, delimiter="\n", skipinitialspace=True)
        data = list(reader)
    ch1 = []
    y1 = []
    aux = []
    for i in range(128):
        aux.append(data[i][0].split())
    #print(aux)
    for i in range(len(aux)):
        for k in range(8):
            ch1.append((int(8*(i)+k+1))*1.8511+69.458) ## axes in keV
            y1.append(float(aux[i][k]))

    ## Creates channel yield lists for File2 and
    # converts channel to energy
    with open(File2, 'r') as file:
        reader = csv.reader(file, delimiter="\n", skipinitialspace=True)
        data = list(reader)
    ch2 = []
    y2 = []
    aux = []
    for i in range(128):
        aux.append(data[i][0].split())
    #print(aux)
    for i in range(len(aux)):
        for k in range(8):
            ch2.append((int(8*(i)+k+1))*1.8511+69.458) ## axes in keV
            y2.append(float(aux[i][k]))
    
    fig, ax = plt.subplots()
    ax.plot(ch1,y1,'.', color ='xkcd:blue', label=(str(File1)))
    ax.plot(ch2,y2,'+', color ='xkcd:red', label=(str(File2)))
    #ax.semilogy(ch1,y1,'.', color ='xkcd:blue', label=(str(File1)))
    #ax.semilogy(ch2,y2,'+', color ='xkcd:red', label=(str(File2)))
    legend = ax.legend(loc="upper right",ncol=1, shadow=False,fancybox=True,framealpha = 0.0,fontsize=20)
    legend.get_frame().set_facecolor('#DAEBF2')
    tick_params(axis='both', which='major', labelsize=22)
    xlabel('Energy (keV)',fontsize=22)
    ylabel('Yield', fontsize=22)
    grid()
    show()
#############################################
#############################################

###                 Plots                     ###
#VdG_Plot('20220907/ERDrun01.dat') # V-Ta-Nb    
VdG_Plot('20220907/ERDrun03.dat') # Yb-Ge-Si-O  
#################################################
#VdG_PlotBoth('20220907/RBS1run01.dat', '20220907/RBS1run02.dat')
#VdG_PlotBoth('20220907/RBS2run01.dat', '20220907/RBS2run02.dat')
#VdG_PlotBoth('20220907/ERDrun01.dat', '20220907/ERDrun02.dat')
#################################################
VdG_PlotBoth('20220907/RBS1run04.dat', '20220907/RBS1run05.dat')
VdG_PlotBoth('20220907/RBS2run04.dat', '20220907/RBS2run05.dat')
VdG_PlotBoth('20220907/ERDrun04.dat', '20220907/ERDrun05.dat')
#################################################