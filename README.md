# Simulation de Dynamique Moléculaire classique avec une méthode Monte-Carlo et l'algorithme de Metropolis

Ce dépôt comprend des codes de dynamique moléculaire écrits en langage C ainsi que des petits codes Python permettant l'analyse complémentaire des données issues des simulations. 

Tous les codes sont commentés. 

## :thread: Contexte

Dans le cadre de l'enseignement de **Simulation Atomistique des Matériaux** des master 2 de Physique Numérique et de Chimie Théorique de l'université de Montpellier les codes suivants ont été développés. Le but étant d'avoir un code de dynamique moléculaire classique 2D robuste utilisant la méthode monte-carlo pour pouvoir illustrer les transitions de phase. 

Le potentiel utilisé ici est le potentiel de Lennard-Jones dans un cadre 2d : $u(r) = 4 \epsilon \left[ (\frac{\sigma}{r})^{12} - (\frac{\sigma}{r})^{6} \right]$. Les N particules de notre système étudié se trouvent dans une boîte 2D pour laquelle on a appliqué les conditions aux bords périodiques. Des paramètres comme le nombre de particules ou encore la taille de la boîte peuvent être modifiés directement dans le code. 

Attention on utilise la formule du potentiel de Lennard-Jonnes adimensionnée par défaut dans le code. 

## :wrench: Méthode de Monte-Carlo - détails

Ces simulations utilisent la méthode de Monte-Carlo mais plus précisément encore selon l'agorithme de Metropolis. Les déplacements aléatoires de particules sont ainsi acceptées ou refusées selon le critère de Metropolis. Chaque tentative de déplacement correspond à une tentative de déplacer **une seule** particule. Et un cycle correspond à N tentatives de déplacement, avec N le nombre de particules dans votre système. Ainsi si vous avez par exemple N=100 et que vous effectuer 5 cycles, vous aurez effectué 500 tentatives, acceptées ou refusées. Ainsi chaque particule, au cours d'un cycle, effectue en moyenne une tentative de déplacement. 

## :spider_web: Configurations initiales 

Il y a sur ce dépôt deux dossiers. 

Le premier nommé "*Configuration initiale aléatoire*" correspond au cas de figure où on initialise aléatoirement la position des N particules au sein de notre boîte au début de la simulation. Cela signifie que si la taille de la boite est de 10, chaque particule aura son x et son y dans des intervalles [-L/2 ; L/2] = [-5;5] les positionnant dans la boîte. 

Le deuxième nommé "*Cristal*" correspond au cas ou initialement notre système est dans une configuration cristalline avec un cristal infiniment périodique. Pour ce cas de figure, si vous voulez avoir réellemnt un cristal infiniment périodique il est sage de choisir un nombre N de particules carré, tel que 6x6=36 ou 10x10=100 particules par exemple. Le programme calculera automatiquement comment espacer vos particules pour que le cristal soit infiniment périodique en prenant en compte les images périodiques. 

<img width="197.5" alt="cristal_infini_1" src="https://github.com/user-attachments/assets/6ba9b531-ef50-4a5b-b935-977eb512bed7">

<img width="335" alt="cristal_infini_distances" src="https://github.com/user-attachments/assets/6d0aead2-b25b-4dd0-a3ef-9e220450840f">

## 🎲 Configuration initiale aléatoire 

Pour ce cas de figure le code C vous sortira un fichier txt nommé "*positions_data.txt*" et qui contiendra les positions en x et y des particules à chaque cycle. Ce fichier est compatible avec Ovito.
Vous aurez aussi dans le terminal l'avancement du nombre de cycle au fur et à mesure de l'exécution ainsi qu'à la fin le taux d'acceptation d'affiché qui correspond à : $taux = \frac{\text{nombre tentatives acceptées}}{\text{nombre total de tentatives}}$. 

Il vous est aussi possible de modifier aisément le code afin de récupérer la pression ou l'énergie potentielle au cours de la simulation, pour cela il suffit de prendre exemple sur les codes C concernant le cristal.

Dans les constantes au début du programme il vous est possible de changer le nombre de particules, la taille de la boîte, le rayon de coupure vis à vis du potentiel, le nombre de cycles, la température ou encore le $\delta_{max}$ qui correspond au rayon de la sphère dans laquelle je déplace ma particule lors d'une tentative. 

## 💎 Configuration initiale cristalline

Pour ce cas-ci il y a 3 dossiers qui vous sont disponibles, chacun servant un objectif différent. 

### Energie potentiel + pression

Ces codes permettent l'étude de la pression et de l'énergie potentiel au cours de la simulation. Ils peuvent permettrent par exemple de déterminer au bout de combien de cycles le système arrive dans un état d'équilibre (plateau de l'énergie potentielle). 
Vous pouvez aussi vous en servir pour voir les influences d'un échauffement ou encore d'un refroidissement de votre système. 

### Isothermes

Ces codes permettent de tracer l'évolution de la pression en fonction du volume du système (en faisant varier la taille de la boîte). Cela permet ainsi de tracer des isothermes P(V). 
Pour calculer la pression la méthode du viriel est utilisée. 

### Transition de phase 

Ces derniers codes permettent de tracer l'évolution de l'énergie potentielle en fonction de la densité de votre système. On fait ainsi évoluer la densité du système et pour chaque densité après que le système soit arrivé à l'équilibre on récupère son énergie potentielle. 
Cela permet d'observer une transition de phase pour par exemple la température kT=0.05. (Avec k la constante de Boltzmann addimensionnée.)

## :ring_buoy: Comment faire fonctionner les codes ? 

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
