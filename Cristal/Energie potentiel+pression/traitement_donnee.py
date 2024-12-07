import matplotlib.pyplot as plt 


fichier="energie_data_alter.txt"
energie=[]
cycles=[]

with open(fichier, 'r') as file : 
    for line in file: 
        parts = line.split(':')

        energie.append(float(parts[1].strip()))
        cycles.append(float(parts[0].strip()))

plt.plot(cycles,energie)
plt.show()