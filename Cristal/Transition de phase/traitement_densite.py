import matplotlib.pyplot as plt 


fichier="energie_densite.txt"
energie=[]
densite=[]

with open(fichier, 'r') as file : 
    for line in file: 
        parts = line.split(':')

        energie.append(float(parts[1].strip()))
        densite.append(float(parts[0].strip()))

plt.plot(densite[-74:],energie[-74:],"o", color="purple")
plt.xlabel("rho (densité)")
plt.ylabel("Energie moyenne (par particule)")
plt.axhline(y=0.0, color='r', linestyle='-')
plt.grid()
plt.show()
