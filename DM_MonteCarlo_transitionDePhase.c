#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>

//-----------------------------{Définition des constantes}-----------------------------
const int nbx_particules = 100;   //Nombre de particules du système
#define Rc 2.5* sigma              //Rayon de coupure 
#define pi M_PI                     //Constante pi
#define L  40               // Taille de la boîte
const int sigma = 1;                // sigma du lennard jones adimentionné
const double Lmoitie = L * 0.5;    // calcul de la valeur de la moitié de L
const double rho = nbx_particules / (L * L); // rho la densité volumique
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
const double U_tail = pi * epsilon * nbx_particules * rho * (terme_commun_frac * inv_Rc_10 - inv_Rc_4);
const double U_decal  = 4 * (inv_Rc_12 - inv_Rc_6);
const double P_tail = 6*pi*epsilon*rho*rho*(terme_commun_frac*inv_Rc_10 - 0.5*inv_Rc_4);

//Je définis le delta max pour le placement aléatoire des particules
const double delta_max = 0.5*sigma; 
const double temperature = 0.05; // imposé dans l'exo kT = 0.01
const int nombre_cycles =1000; //1 cycle = nbx_particules tentatives, réussies ou non

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


void initialiser_Configuration_Cristalline(Particule tab_par[]){

    //Calculer le nombre de particules par ligne et par colonne
    int particules_par_ligne = (int) sqrt(nbx_particules); 
    int particules_par_colonne = (nbx_particules+particules_par_ligne-1)/particules_par_ligne; 

    //Espacement entre les particules
    double espacement_x = 4.0; 
    double espacement_y = 4.0; 

    //Calcul de la largeur et hauteur du réseau 
    double largeur_reseau = espacement_x*(particules_par_ligne-1); 
    double hauteur_reseau = espacement_y*(particules_par_colonne-1); 

    //Décalage pour centrer le réseau atour de (0,0) 
    double decalage_x = -largeur_reseau/2.0; 
    double decalage_y = -hauteur_reseau/2.0; 

    for (int i=0; i<= particules_par_ligne; i++){
        for (int j=0; j< particules_par_colonne; j++){

            int index = i*particules_par_colonne+j; 
            if (index<nbx_particules){

                // Position initiale centrée autour de (0, 0)
                tab_par[index].x = i * espacement_x + decalage_x;
                tab_par[index].y = j * espacement_y + decalage_y;

                tab_par[index].actif = 1;
            }
        }
    }
}


//-----------------------------{changement de configuration aleatoire}------------------------------------------
/** 
 * 
 * 
*/
void change_nouvelle_configuration(Particule tab_par[], int nbx_particules)
{
    // Générer un nombre entier aléatoire entre 0 et N exclus
    int index_particule = rand() % (nbx_particules);

    // Générer un rayon aléatoire entre 0 et delta_max
    double rayon = ((double)rand() / RAND_MAX) * delta_max;

    // Générer un angle aléatoire entre 0 et 2*pi (en radians)
    double angle = ((double)rand() / RAND_MAX) * 2 * M_PI;

    // Calculer les déplacements en x et y en utilisant le rayon et l'angle
    double delta_x = rayon * cos(angle);
    double delta_y = rayon * sin(angle);

    // Ajouter ces déplacements à la position actuelle de la particule
    tab_par[index_particule].x += delta_x;
    tab_par[index_particule].y += delta_y;

    //Application des conditions périodiques sur x : 
    if (tab_par[index_particule].x >= Lmoitie) {tab_par[index_particule].x -=L;}
    if (tab_par[index_particule].x < -Lmoitie) {tab_par[index_particule].x += L;}

    //Application des conditions périodiques sur y : 
    if (tab_par[index_particule].y >= Lmoitie) {tab_par[index_particule].y -=L;}
    if (tab_par[index_particule].y < -Lmoitie) {tab_par[index_particule].y += L;}

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

    if (energie_tab_temp <= energie_tab_ini){
        proba = 1.0; 
    } else {
        proba = exp(-(energie_tab_temp-energie_tab_ini)/(temperature));
    }

    return proba; 
}



//-----------------------------{Enregistrer les positions}-------------------------------------------

void enregistrer_positions(FILE *fp, Particule *tab_par, int tentative, int nbx_particules){
    fprintf(fp, "%d\n", nbx_particules);
    fprintf(fp, "\n");
    for (int i=0; i<nbx_particules; i++){
        fprintf(fp, "A %f %f\n", tab_par[i].x, tab_par[i].y);
    }
}

void enregistrer_energie(FILE *fp, double energie, int cycle){
    fprintf(fp, "%d : %f\n", cycle, energie); 
}

//-----------------------------{Vérifier les différences}-------------------------------------------

/**
 * Fonction qui permet de vérifier quelle particule a bougée entre deux configurations
 */
void verifier_differences(Particule tab_ini[], Particule tab_temp[], int nbx_particules) {
    for (int i = 0; i < nbx_particules; i++) {
        if (tab_ini[i].x != tab_temp[i].x || tab_ini[i].y != tab_temp[i].y) {
            printf("Particule %d a changé : \n", i);
            printf("  Initial : x = %f, y = %f\n", tab_ini[i].x, tab_ini[i].y);
            printf("  Temporaire : x = %f, y = %f\n", tab_temp[i].x, tab_temp[i].y);
        }
    }
}



//-----------------------------{main}----------------------------------------------------------------
int main() 
{
    // Initialiser le générateur de nombres aléatoires avec une graine pour avoir des résultats reproductibles
    srand(42);

    // Creation d'un tableau de particules
    Particule tab_particules[nbx_particules];
    Particule tab_particules_temporaire[nbx_particules]; 

    double energie_totale[nombre_cycles+1]; // Je cherche à voir au bout de combien de cycles j'atteinds l'équilibre

    // Initialisation des particules avec une configuration aleatoire
    initialiser_Configuration_Cristalline(tab_particules); 
    
    // Ouvre le fichier "positions_data.txt" en mode écriture ("w")
    // Ce fichier stockera les positions des particules au fil du temps dans le format suivant:
    // temps   x1 y1    x2 y2    x3 y3    x4 y4.....
    FILE *pos_file = fopen("positions_data.txt", "w");
    FILE *energie_file = fopen("energie_data.txt", "w"); 

    int tentative = 0;
    int tentative_acceptee=0; 

    // Appelle la fonction enregistrer_Positions pour écrire les positions initiales des particules
    enregistrer_positions(pos_file, tab_particules, tentative, nbx_particules);

    for (int cycle=0 ; cycle<nombre_cycles; cycle++){

            double U = calcul_energie_potentielle(tab_particules); 
            energie_totale[cycle] = U; 
            enregistrer_energie(energie_file, U, cycle); 

            for (int tentativou=0; tentativou < nbx_particules ; tentativou++){

                for (int i=0; i<nbx_particules; i++){
                    tab_particules_temporaire[i]=tab_particules[i];
                }

                change_nouvelle_configuration(tab_particules_temporaire, nbx_particules); 

                // Vérifier les différences entre les deux configurations
                //verifier_differences(tab_particules, tab_particules_temporaire, nbx_particules);


                double proba_finale=calcul_boltzmann(tab_particules, tab_particules_temporaire); 
                //printf("%f\n", proba_finale);

                double x= ((double)rand() / RAND_MAX);

                if (proba_finale >= x){

                    tentative_acceptee += 1; 

                    for (int i=0; i<nbx_particules; i++){
                        tab_particules[i]=tab_particules_temporaire[i];
                    }
                }

                tentative += 1; 
                //printf("%d\g", tentative_acceptee);
            }

            enregistrer_positions(pos_file, tab_particules, tentative, nbx_particules); 
    }

    double U = calcul_energie_potentielle(tab_particules); 
    energie_totale[nombre_cycles]; 
    enregistrer_energie(energie_file, U, nombre_cycles); 

    fclose(pos_file);  // De même pour le fichier positions_data.txt

    int nbx_total_tentative = nbx_particules*nombre_cycles; 
    double taux_acceptation = (double)tentative_acceptee / nbx_total_tentative;
    printf("Taux d'acceptation : %g\n", taux_acceptation); 

    return 0;
}
