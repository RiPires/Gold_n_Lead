#RiP
from matplotlib.pylab import *
import matplotlib.pyplot as plt
import csv
from scipy import stats



##   Plots spectrum from DATA file   ###############################
def PlotData(File):
    """
    Plots yield vs energy data from .csv alfaMC_3body output files

    INPUTS: "FileName.csv" output from alfaMC
    OUTPUTS: 
    """

    with open(File, 'r') as file:
        reader = csv.reader(file, delimiter="\n", skipinitialspace=True)
        data = list(reader)
    
    x = []
    y = []
    aux = []
    
    for i in range(2,len(data)): ## 2 is the third line
        aux.append(data[i][0].split())
    for i in range(len(aux)):
        x.append(float(aux[i][0]))
        y.append(float(aux[i][1]))

    fig, ax = plt.subplots()
    ax.plot(x,y,'-', color ='xkcd:black', label=(str(File)))
    #ax.semilogy(x,y,'^', color ='xkcd:purple', label=(str(File)))
    legend = ax.legend(loc="upper right",ncol=1, shadow=True,fancybox=True,framealpha = 0.0,fontsize=20)
    legend.get_frame().set_facecolor('#DAEBF2')
    tick_params(axis='both', which='major', labelsize=22)
    xlabel('Energy',fontsize=22)
    ylabel('Counts', fontsize=22)
    show()

    return
########################################################################

##   Plots spectrum from SIMULATION file   #############################
def PlotSimu(File):

    with open(File, 'r') as file:
        reader = csv.reader(file, delimiter="\n", skipinitialspace=True)
        data = list(reader)
    
    x = []
    y = []
    aux = []
    
    for i in range(2,len(data)): ## 2 is the third line
        aux.append(data[i][0].split())
    for i in range(len(aux)):
        x.append(float(aux[i][0]))
        y.append(float(aux[i][1]))

    fig, ax = plt.subplots()
    ax.plot(x,y,'-', color ='xkcd:black', label=(str(File)))
    #ax.semilogy(x,y,'^', color ='xkcd:black', label=(str(File)))
    legend = ax.legend(loc="upper right",ncol=1, shadow=True,fancybox=True,framealpha = 0.0,fontsize=20)
    legend.get_frame().set_facecolor('#DAEBF2')
    tick_params(axis='both', which='major', labelsize=22)
    xlabel('Energy',fontsize=22)
    ylabel('Counts', fontsize=22)
    show()

    return
########################################################################

##   Plots both spectrums from data and simulation   ###############################
def PlotBoth(dataFile, simuFile):

    ### opens data file and fills lists   ###
    with open(dataFile, 'r') as fileData:
        reader = csv.reader(fileData, delimiter="\n", skipinitialspace=True)
        data = list(reader)
    data_Energy = []
    data_Counts = []
    aux = []
    for i in range(12,len(data)-1): ## 12
        aux.append(data[i][0].split())
    for i in range(len(aux)):
        data_Energy.append(float(0.004482*i-0.105144)) ## convert channel to energy, calibration updated
        data_Counts.append(float(aux[i][0]))

    ###   opens simulation file and fills lists   ###
    with open(simuFile, 'r') as fileSimu:
        reader = csv.reader(fileSimu, delimiter="\n", skipinitialspace=True)
        dataSimu = list(reader)
    simu_Energy = [] ## in keV
    simu_Counts = []
    aux = []
    for i in range(2,len(dataSimu)): ## 2 is the third line
        aux.append(dataSimu[i][0].split())
    for i in range(len(aux)):
        simu_Energy.append(float(aux[i][0]))
        simu_Counts.append(float(aux[i][1]))

    max_data = max(data_Counts)
    max_simu = max(simu_Counts)

    for k in range(len(data_Counts)):
        data_Counts[k] = data_Counts[k]/max_data
    for k in range(len(simu_Counts)):
        simu_Counts[k] = simu_Counts[k]/max_simu
    

    ###   Plots   ###
    fig, ax = plt.subplots()
    ax.plot(data_Energy,data_Counts,'.', color ='xkcd:black', label=(str(dataFile)))
    ax.plot(simu_Energy,simu_Counts,'-', color ='xkcd:red', label=(str(simuFile)))
    legend = ax.legend(loc="upper left",ncol=1, shadow=True,fancybox=True,framealpha = 0.0,fontsize=20)
    legend.get_frame().set_facecolor('#DAEBF2')
    tick_params(axis='both', which='major', labelsize=22)
    xlabel('Energy [keV]',fontsize=22)
    ylabel('Normalized yield', fontsize=22)
    show()

    return
########################################################################


##   Plots both spectrums from different data   ###############################
def PlotBothData(dataFile1, dataFile2):

    ### opens data file 1 and fills lists   ###
    with open(dataFile1, 'r') as fileData1:
        reader = csv.reader(fileData1, delimiter="\n", skipinitialspace=True)
        data = list(reader)
    data_Energy1 = []
    data_Counts1 = []
    aux = []
    for i in range(12,len(data)-1): ## 2 is the third line
        aux.append(data[i][0].split())
    for i in range(len(aux)):
        data_Energy1.append(float(4.5533*i-3.8118)) ## convert channel to energy, calibration updated
        data_Counts1.append(float(aux[i][0]))

    ### opens data file 1 and fills lists   ###
    with open(dataFile2, 'r') as fileData2:
        reader = csv.reader(fileData2, delimiter="\n", skipinitialspace=True)
        data = list(reader)
    data_Energy2 = []
    data_Counts2 = []
    aux = []
    for i in range(12,len(data)-1): ## 2 is the third line
        aux.append(data[i][0].split())
    for i in range(len(aux)):
        data_Energy2.append(float(4.5533*i-3.8118)) ## convert channel to energy, calibration updated
        data_Counts2.append(float(aux[i][0]))

    max_data1 = max(data_Counts1)
    max_data2 = max(data_Counts2)

    for k in range(len(data_Counts1)):
        data_Counts1[k] = data_Counts1[k]/max_data1
    for k in range(len(data_Counts2)):
        data_Counts2[k] = data_Counts2[k]/max_data2
    
    ###   Plots   ###
    fig, ax = plt.subplots()
    ax.plot(data_Energy1,data_Counts1,'.', color ='xkcd:black', label=(str(dataFile1)))
    ax.plot(data_Energy2,data_Counts2,'-', color ='xkcd:red', label=(str(dataFile2)))
    legend = ax.legend(loc="upper left",ncol=1, shadow=True,fancybox=True,framealpha = 0.0,fontsize=20)
    legend.get_frame().set_facecolor('#DAEBF2')
    tick_params(axis='both', which='major', labelsize=22)
    xlabel('Energy [keV]',fontsize=22)
    ylabel('Normalized yield', fontsize=22)
    show()

    return
########################################################################

###   Plot stopping power from NIST   ###
def PlotStopPow(File):
    """
    File .csv in two column: E StopPow
    """
    with open(File, 'r') as file:
        reader = csv.reader(file, delimiter="\n", skipinitialspace=True)
        data = list(reader)
    
    Energy = []
    StopPow = []
    aux = []
    
    for i in range(0,len(data)): ## 0 is the first line
        aux.append(data[i][0].split())
    for i in range(len(aux)):
        Energy.append(float(aux[i][0]))
        StopPow.append(float(aux[i][1]))

    ###   Lower energy region linear extrapolation   ###
    StopLin = []
    EnergyLin = []
    for i in range(int((4.75-4.5)/0.01)):
        StopLin.append(StopPow[i])
        EnergyLin.append(Energy[i])

    slope, intercept, r_value, p_value, std_err = stats.linregress(EnergyLin, StopLin)   

    StopLin = []
    EnergyLin = []
    for k in range(len(Energy)):
        EnergyLin.append(Energy[k])
        StopLin.append(slope * Energy[k] + intercept)
    ###################################################################

    ###   Higher energy region linear extrapolation   ###
    StopLin2 = []
    EnergyLin2 = []
    for i in range(int(7.25*320/7.7), 320, 1):
        StopLin2.append(StopPow[i])
        EnergyLin2.append(Energy[i])

    slope2, intercept2, r_value, p_value, std_err = stats.linregress(EnergyLin2, StopLin2)   

    StopLin2 = []
    EnergyLin2 = []
    for k in range(len(Energy)):
        EnergyLin2.append(Energy[k])
        StopLin2.append(slope2 * Energy[k] + intercept2)
    ###################################################################

    fig, ax = plt.subplots()
    ax.plot(Energy,StopPow,'-', color ='xkcd:blue', label=(str(File)))
    #ax.plot(EnergyLin,StopLin,'--', color ='xkcd:red')
    ax.plot(EnergyLin2,StopLin2,'--', color ='xkcd:purple')
    legend = ax.legend(loc="upper right",ncol=1, shadow=False,fancybox=True,framealpha = 0.0,fontsize=20)
    legend.get_frame().set_facecolor('#DAEBF2')
    tick_params(axis='both', which='major', labelsize=22)
    xlabel('Energy (MeV)',fontsize=22)
    ylabel('Stop. Pow. (Mev cm$^2$ g$^{-1}$)', fontsize=22)
    show()

    return
    
########################################################################

###   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   ###


PlotBoth('u232+Au.mca', 'alfaMC_filme_Au/Edet.csv')