import csv

Output = open('U232_Calibrated_v1.in', 'w')

with open('U232_v1.in', 'r') as file:
    reader = csv.reader(file, delimiter="\n", skipinitialspace=True)
    data = list(reader) 
counts = []
aux = []

for i in range(0,len(data)):
    aux.append(data[i][0].split())
for i in range(len(aux)):
    counts.append(int(aux[i][0]))

for i in range(len(counts)):
    #Output.write("{:.2f}".format((i+1)*0.004941-0.10356)+'\t'+str(counts[i])+'\n')
    #Output.write("{:.2f}".format((i+1)*0.00501-0.025)+'\t'+str(counts[i])+'\n') ## (not so) better calib. parameters
    Output.write("{:.2f}".format((0.004482*(i+1)-0.105144))+'\t'+str(counts[i])+'\n')

Output.close()