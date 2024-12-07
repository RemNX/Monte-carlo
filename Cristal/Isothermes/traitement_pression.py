import matplotlib.pyplot as plt 


fichier="pression_data_1.txt"
pression_1=[]
volume=[]

with open(fichier, 'r') as file : 
    for line in file: 
        parts = line.split(':')

        pression_1.append(float(parts[1].strip()))
        volume.append(float(parts[0].strip()))

pression_075=[]
pression_05=[]
pression_025=[]
pression_125=[]
pression_15=[]

fichier2="pression_data_0-75.txt"
fichier3="pression_data_0-5.txt"
fichier4="pression_data_0-25.txt"
fichier5="pression_data_1-25.txt"
fichier6="pression_data_1-5.txt"

with open(fichier2, 'r') as file : 
    for line in file: 
        parts = line.split(':')
        pression_075.append(float(parts[1].strip()))

with open(fichier3, 'r') as file : 
    for line in file: 
        parts = line.split(':')
        pression_05.append(float(parts[1].strip()))

with open(fichier4, 'r') as file : 
    for line in file: 
        parts = line.split(':')
        pression_025.append(float(parts[1].strip()))

with open(fichier5, 'r') as file : 
    for line in file: 
        parts = line.split(':')
        pression_125.append(float(parts[1].strip()))

with open(fichier6, 'r') as file : 
    for line in file: 
        parts = line.split(':')
        pression_15.append(float(parts[1].strip()))


plt.plot(volume,pression_1,label="kT=1")
plt.plot(volume,pression_075,label="kT=0.75")
plt.plot(volume,pression_05,label="kT=0.5")
plt.plot(volume,pression_025,label="kT=0.25")
plt.plot(volume,pression_125,label="kT=1.25")
plt.plot(volume,pression_15,label="kT=1.5")
plt.xlabel("Volume (m²)")
plt.ylabel("P (S.I.)")
plt.legend()
plt.show()
