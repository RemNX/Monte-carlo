# Simulation de Dynamique Moléculaire classique avec une méthode Monte-Carlo et l'algorithme de Metropolis

Ce dépôt comprend des codes de dynamique moléculaire écrits en langage C ainsi que des petits codes Python permettant l'analyse complémentaire des données issues des simulations. 

## Contexte

Dans le cadre de l'enseignement de **Simulation Atomistique des Matériaux** des master 2 de Physique Numérique et de Chimie Théorique de l'université de Montpellier les codes suivants ont été développés. Le but étant d'avoir un code de dynamique moléculaire classique 2D robuste utilisant la méthode monte-carlo pour pouvoir illustrer les transitions de phase. 

Le potentiel utilisé ici est le potentiel de Lennard-Jones dans un cadre 2d : $u(r) = 4 \epsilon \left[ (\frac{\sigma}{r})^{12} - (\frac{\sigma}{r})^{6} \right]$. Les N particules de notre système étudié se trouvent dans une boîte 2D pour laquelle on a appliqué les conditions aux bords périodiques. Des paramètres comme le nombre de particules ou encore la taille de la boîte peuvent être modifiés directement dans le code. 


## Comment faire fonctionner les codes ? 

Afin de faire fonctionner les codes se trouvant sur le github vous devez suivre les étapes suivantes : 
- Choisir si vous souhaitez étudier le cas d'une configuration initiale aléatoire ou cristalline.
- Dans le cas de la configuration initiale cristalline choisir quelles sont vos grandeurs d'intérêts (pression en fonction du volume, énergie potentielle au cours de la simulation, etc...).
- Une fois vos choix faits grâce à ce petit README vous devriez être en mesure de déterminer quel est le code C qui vous intéresse. Télécharger le fichier C sur votre ordinateur.
- Ouvrez ensuite votre terminal. Si vous n'êtes pas sous un système d'exploitation UNIX vous pouvez aussi ouvrir WSL dans notre système d'exploitation Windows par exemple.
- Naviguez à l'aide des commandes bash **cd** et **ls** jusqu'au répertoire où vous avez téléchargé le code C. 
- Compilez ensuite le code C avec le flag -lm. Vous pouvez aussi mettre des flags d'optimisation comme par exemple -O2 ou -O3, ils ont été testés et n'influencent pas la validité des résultats.
  Par exemple sur un ordinateur sous Windows dans WSL je peux compiler de cette manière :
  **gcc -O2 -o main_DM_simulation main_DM_simulation.c -lm**
- Exécutez ensuite le code C.
  Pour reprendre notre exemple ce sera alors : **time ./main_DM_Simulation**
- Une fois le code exécuté un fichier txt contenant vos données devrait normalement apparaître dans votre répertoire. Vous pouvez alors utilisé le fichier python (préalablement téléchargé dans le même répertoire) correspondant au code C pour analyser vos données (obtenir des graphes surtout).
- Si le fichier de données obtenu avec votre code C contient les positions des particules il est normalement sous le bon format pour être utilisé sur Ovito et vous pourrez alors voir l'évolution du système en fonction du nombre de cycles.  
