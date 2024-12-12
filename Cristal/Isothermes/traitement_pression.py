import matplotlib.pyplot as plt 


fichier="pression_data_5.txt"
pression_5=[]
volume=[]

with open(fichier, 'r') as file : 
    for line in file: 
        parts = line.split(':')
        pression_5.append(float(parts[1].strip()))
        volume.append(float(parts[0].strip()))


plt.plot(volume,pression_5,"o",label="kT=1")
plt.xlabel("Volume (m²)")
plt.ylabel("P (S.I.)")
plt.legend()
plt.show()
