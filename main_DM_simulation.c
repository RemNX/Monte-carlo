#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>

//-----------------------------{Définition des constantes}-----------------------------
const int nbx_particules = 128;   //Nombre de particules du système
const int nbx_particules_actives = 128; //Nombre de particules actives du système (sert juste pour définir des défauts)
#define Rc 2.5* sigma              //Rayon de coupure 
#define pi M_PI                     //Constante pi
#define L  11               // Taille de la boîte
const int sigma = 1;                // sigma du lennard jones adimentionné
const double Lmoitie = L * 0.5;    // calcul de la valeur de la moitié de L
const double rho = nbx_particules_actives / (L * L); // rho la densité volumique
const double epsilon = 1;          //le E0 du lennard jones adimentionné

/*Calcul pour Rc pour eviter de les refaire, vu qu'ils sont constants*/
const double Rc2 = Rc * Rc;        
const double inv_Rc_2 = 1 / Rc2;
const double inv_Rc_4 = inv_Rc_2 * inv_Rc_2;
const double inv_Rc_6 = inv_Rc_2 * inv_Rc_2 * inv_Rc_2;
const double inv_Rc_10 = inv_Rc_6 * inv_Rc_4;
const double inv_Rc_12 = inv_Rc_6 * inv_Rc_6;

/*Calcul des terme à ajouter pour les grandeurs calculées, qui sont constants*/
const double terme_commun_frac = 2/5;
const double U_tail = pi * epsilon * nbx_particules_actives * rho * (terme_commun_frac * inv_Rc_10 - inv_Rc_4);
const double U_decal  = 4 * (inv_Rc_12 - inv_Rc_6);
const double P_tail = 6*pi*epsilon*rho*rho*(terme_commun_frac*inv_Rc_10 - 0.5*inv_Rc_4);

//Je définis le delta max pour le placement aléatoire des particules
const double delta_max = 0.2*sigma; 
const double temperature = 0.01; // imposé dans l'exo kT = 0.01
const int nombre_cycles =10; //1 cycle = 128 tentatives, réussies ou non

//-----------------------------{Définition de la structure pour la particule}-----------------------------
/**
 * Définition de la classe Particule avec les attribus.
 *
 * @param x Coordonnée x de la particule.
 * @param y Coordonnée y de la particule
 * @param vx vitesse selon x de la particule.
 * @param vy vitesse selon y de la particule.
 * @param fx force appliquée selon x sur la particule.
 * @param fy force appliquée selon y sur la particule.
 * @param actif  indique si la particule existe(sert pour les défauts)
 */
typedef struct {
    double x; 
    double y;
    int actif; //1=active et 0= inactive, permet d'inclure des défauts
} Particule;

//-----------------------------{Constructeur pour la particule}---------------------------------------
/**
 * @brief Initialise les attributs d'une instance de Particule.
 *
 * Cette fonction configure les coordonnées et les vitesses d'une particule 
 * en initialisant ses attributs avec les valeurs fournies.
 *
 * @param p Pointeur vers l'instance de Particule à initialiser.
 * @param x Coordonnée x de la particule.
 * @param y Coordonnée y de la particule.
 */
void init_Particule(Particule *p, double x, double y) {
    // Initialisation des coordonnées de la particule
    p->x = x;
    p->y = y;
}

/**
 * @brief Initialise une configuration aléatoire de particules dans une boîte carrée de côté L.
 * 
 * Cette fonction remplit un tableau de particules avec des positions et des vitesses générées aléatoirement.
 * Les particules sont positionnées dans une plage de coordonnées allant de -L/2 à L/2.
 * Chaque nouvelle particule respecte une distance minimale spécifiée par rapport aux particules existantes.
 * 
 * @param tab_par Tableau de particules à remplir avec des configurations aléatoires.
 */
void initialiser_Configuration_aleatoire(Particule tab_par[]) {
    int fac = 0; // Facteur maximum pour les valeurs aléatoires de vitesse (vx, vy)
    double MIN_DISTANCE = 1; // Distance minimale entre les particules pour éviter les chevauchements

    // Boucle pour générer nbx_particules particules
    for (int i = 0; i < nbx_particules; i++) {
        Particule p; // Nouvelle particule à ajouter au tableau
        int isTooClose; // Flag pour vérifier si la particule est trop proche d'une autre

        // Génère une nouvelle particule jusqu'à ce qu'elle soit suffisamment éloignée des autres
        do {
            isTooClose = 0; // Réinitialise le flag avant de générer une nouvelle particule

            // Génère des coordonnées (x, y) aléatoires entre -L/2 et L/2
            p.x = ((double)rand() / RAND_MAX * L) - Lmoitie;
            p.y = ((double)rand() / RAND_MAX * L) - Lmoitie;

            // Vérifie la distance de cette particule avec chaque particule déjà générée
            for (int j = 0; j < i; j++) {
                // Calcul de la distance entre la particule actuelle et chaque particule déjà existante
                double distance = sqrt(pow(p.x - tab_par[j].x, 2) + pow(p.y - tab_par[j].y, 2));
                
                // Si la particule est trop proche d'une autre, on active le flag et on quitte la boucle
                if (distance < MIN_DISTANCE) {
                    isTooClose = 1; // Marque que la particule est trop proche
                    break; // Sort de la boucle de vérification
                }
            }
        } while (isTooClose); // Répète jusqu'à trouver une position non conflictuelle

        // Ajoute la particule validée dans le tableau
        tab_par[i] = p;
        tab_par[i].actif=1;
    }
}



//-----------------------------{changement de configuration aleatoire}------------------------------------------
/** 
 * 
 * 
*/
void change_nouvelle_configuration(Particule tab_par[], int nbx_particules)
{
    //Je choisis l'index aléatoirement pour choisir une particule aléatoirement dans ma configuration
    srand(time(NULL)); //! A améliorer ensuite
    // Générer un nombre entier aléatoire entre 0 et N exclus
    int index_particule = rand() % (nbx_particules);

    //Maintenant je veux fixer de nouvelles positions aléatoirement :

    double delta_moitie = delta_max/2;

    bool ok_cercle=false;  

    while (ok_cercle == false) {

        tab_par[index_particule].x = ((double)rand() / RAND_MAX *2* delta_max) - delta_max; //position générée entre -L/2 et L/2
        tab_par[index_particule].y = ((double)rand() / RAND_MAX *2* delta_max) - delta_max;

        double x_double = tab_par[index_particule].x*tab_par[index_particule].x; 
        double y_double = tab_par[index_particule].y*tab_par[index_particule].y; 

        if ((x_double + y_double) < delta_max*delta_max){
            ok_cercle=true;
        }

    } 

}

//-----------------------------{Calcul de Ep}-----------------------------------------
/**
 * @brief Calcule les grandeurs thermodynamiques du système.
 *
 * Cette fonction effectue les calculs pour déterminer la température, 
 * la pression, l'énergie potentielle et l'énergie totale du système.
 *
 * @param tab_par Tableau des particules du système.
 * @param temperature Pointeur vers la variable où sera stockée la température calculée.
 * @param energie_potentielle Pointeur vers la variable où sera stockée l'énergie potentielle calculée.
 */
double calcul_energie_potentielle(Particule tab_par[]) 
{
    /*définition des variables locales à utiliser dans les calculs*/
    double U = 0.0; // energie potentielle
    double rij2, inv_r2, inv_r2x3, inv_r2x6; //puissance de rij
    double V = L * L; //volume
    int dimension = 2;

    for (int i = 0; i < nbx_particules; i++) 
    {
        /*condition qui verifie que le calcul se fait que par rapports aux partiucules et ignore les défauts*/
        if (tab_par[i].actif) 
        {
            for (int j = i + 1; j < nbx_particules; j++) 
            {
                if (!tab_par[j].actif) continue; 

                double x_sous = tab_par[j].x - tab_par[i].x; //Dx
                double y_sous = tab_par[j].y - tab_par[i].y; //
                /*Calcul de la distance entre deux particules*/
                rij2 = x_sous * x_sous + y_sous * y_sous;   

                /*Calcul sur les puissances de rij pour eviter l'appel de pow*/
                inv_r2 = 1.0 / rij2;                  //   1/rij^2
                inv_r2x3 = inv_r2 * inv_r2 * inv_r2; //    1/rji^6
                inv_r2x6 = inv_r2x3 * inv_r2x3;     //    1/rij^12

                /*Conditions qui permet de tronquer le potentiel à Rc*/
                if (rij2 < Rc2) 
                {
                    U += 4 * (inv_r2x6 - inv_r2x3) - U_decal; //calcul du potentiel entre deux particules
                }
            }
        }
        
    }
    return U+U_tail;    //calcul et mise à jour l'energie potentielle
}

//-----------------------------{calcul de boltzmann}------------------------------------------

double calcul_boltzmann(Particule tab_ini[], Particule tab_temp[])
{
    double proba = 0.0;

    double energie_tab_ini = calcul_energie_potentielle(tab_ini);
    double energie_tab_temp = calcul_energie_potentielle(tab_temp);

    if (energie_tab_temp >= energie_tab_ini){
        proba = 1.0; 
    } else {
        proba = exp(-(energie_tab_temp-energie_tab_ini)/(temperature));
    }

}


// //-----------------------------{Initialise la configuration cristalline}------------------------------------------
// /**
//  * @brief Initialise une configuration cristalline pour un ensemble de particules.
//  *
//  * Cette fonction dispose les particules dans une structure cristalline,
//  * typiquement en utilisant une grille régulière dans l'espace. 
//  * Cela signifie que chaque particule est placée à une position fixe et ordonnée.
//  *
//  * @param tab_par Pointeur vers un tableau de structures `Particule`, 
//  *                représentant les particules du système à initialiser.
//  *                Le tableau doit être préalablement alloué avec le nombre
//  *                de particules souhaité.
//  *
//  * @note Cette fonction suppose que le tableau `tab_par` est suffisamment 
//  *       grand pour contenir toutes les particules nécessaires.
//  * 
//  */
// void initialiser_Configuration_Cristalline(Particule *tab_par) 
// {
//     int fac = 10; // Facteur maximum pour les valeurs aléatoires de vitesse (vx, vy)
//     // Calculer le nombre de particules par ligne et par colonne
//     int particules_par_ligne = (int)sqrt(nbx_particules); 
//     int particules_par_colonne = (nbx_particules + particules_par_ligne - 1) / particules_par_ligne; 
    
//     // Espacement entre les particules
//     double espacement_x = 1.1; // Ajustez pour obtenir l'espacement souhaité
//     double espacement_y = 1.1;

//     // Calcul de la largeur et hauteur du réseau
//     double largeur_reseau = espacement_x * (particules_par_ligne - 1);
//     double hauteur_reseau = espacement_y * (particules_par_colonne - 1);

//     // Décalage pour centrer le réseau autour de (0, 0)
//     double decalage_x = -largeur_reseau / 2.0;
//     double decalage_y = -hauteur_reseau / 2.0;

//     for (int i = 0; i <= particules_par_ligne; i++) 
//     {
//         for (int j = 0; j < particules_par_colonne; j++) 
//         {
//             int index = i * particules_par_colonne + j;
//             if (index < nbx_particules) 
//             {
//                 // Position initiale centrée autour de (0, 0)
//                 tab_par[index].x = i * espacement_x + decalage_x;
//                 tab_par[index].y = j * espacement_y + decalage_y;
//                 //tab_par[index].vx = 0; 
//                 //tab_par[index].vy = 0;
//                 tab_par[index].fx=0;
//                 tab_par[index].fy=0;

//                 // Générer une vitesse initiale entre -0.5 et 0.5, puis la multiplier par le facteur
//                 tab_par[index].vx = ((double)rand() / RAND_MAX - 0.5) * fac;
//                 tab_par[index].vy = ((double)rand() / RAND_MAX - 0.5) * fac;

//                 tab_par[index].actif = 1; 
//                 //si jamais on veut ajouter des defauts 
//                 /*if(index==12 || index==85)
//                 {
//                     tab_par[index].actif=0; 
//                 }*/
//             }
//         }
//     }
// }


//-----------------------------{Enregistrer les positions}-------------------------------------------
void enregistrer_Positions(Particule *tab_par, double tentative, FILE *file) {
    fprintf (file, "%22.8d",tentative);
        for ( int i=0;i<nbx_particules;i++)
        {

        if (!tab_par[i].actif) continue;

         fprintf (file, "%22.8g %22.8g",tab_par[i].x,tab_par[i].y);
        }
        fprintf (file, "\n");
}


//-----------------------------{main}----------------------------------------------------------------
int main() 
{
    // Initialiser le générateur de nombres aléatoires avec une graine pour avoir des résultats reproductibles
    srand(42);

    // Creation d'un tableau de particules
    Particule tab_particules[nbx_particules];
    Particule tab_particules_temporaire[nbx_particules]; 

    // Initialisation des particules avec une configuration aleatoire
    initialiser_Configuration_aleatoire(tab_particules);
    
    // Ouvre le fichier "positions_data.txt" en mode écriture ("w")
    // Ce fichier stockera les positions des particules au fil du temps dans le format suivant:
    // temps   x1 y1    x2 y2    x3 y3    x4 y4.....
    FILE *pos_file = fopen("positions_data.txt", "w");

    int tentative = 0;

    // Appelle la fonction enregistrer_Positions pour écrire les positions initiales des particules
    enregistrer_Positions(tab_particules, tentative, pos_file);

    for (int cycle=0 ; cycle<nombre_cycles; cycle++){

            for (int tentativou=0; tentativou < nbx_particules ; tentativou++){

                for (int i=0; i<nbx_particules; i++){
                    tab_particules_temporaire[i]=tab_particules[i];
                }

                change_nouvelle_configuration(tab_particules_temporaire, nbx_particules); 

                double proba_finale=calcul_boltzmann(tab_particules, tab_particules_temporaire); 

                double x= ((double)rand() / RAND_MAX);

                if (proba_finale >= x){

                    for (int i=0; i<nbx_particules; i++){
                        tab_particules[i]=tab_particules_temporaire[i];
                    }
                }

                tentative += 1; 

                enregistrer_Positions(tab_particules, tentative, pos_file); 
            }
    }

    fclose(pos_file);  // De même pour le fichier positions_data.txt

    // Ouvre un processus Gnuplot en mode écriture ("w") pour envoyer des commandes de tracé
    // "gnuplot -persistent" permet de garder la fenêtre du graphique ouverte après l'exécution des commandes
    /*FILE *gnuplot = popen("gnuplot -persistent", "w");

    if (gnuplot) {

    // Sous-graphe pour l'énergie potentielle (bleu)
    fprintf(gnuplot, "set ylabel '{/Times-Italic U}' font ',18'\n");
    fprintf(gnuplot, "set title '{/Times-Italic Énergie Potentielle}' font ',20' \n");
    fprintf(gnuplot, "plot 'simulation_data.txt' using 1:3 with lines linecolor rgb '#1976D2' title 'Énergie Potentielle', \\\n");
    fprintf(gnuplot, "     (%.6f) title 'Moyenne' with lines linestyle 2 linecolor rgb '#FF9800'\n", moyenne_ep);

    // Fermeture du processus Gnuplot
    pclose(gnuplot);*/

    return 0;
}
