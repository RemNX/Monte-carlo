#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>

//-----------------------------{Définition des constantes}-----------------------------
const int nbx_particules = 128;   //Nombre de particules du système
const int nbx_particules_actives = 128; //Nombre de particules actives du système (sert juste pour définir des défauts)
const double sigma = 1;                // sigma du lennard jones adimentionné
#define Rc 2.5* sigma              //Rayon de coupure 
#define pi M_PI                     //Constante pi
#define L  40              // Taille de la boîte
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
const int nombre_cycles =100; //1 cycle = 128 tentatives, réussies ou non

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

    // Boucle pour générer nbx_particules particules
    for (int i = 0; i < nbx_particules; i++) {
        Particule p; // Nouvelle particule à ajouter au tableau


        // Génère des coordonnées (x, y) aléatoires entre -L/2 et L/2
        p.x = ((double)rand() / RAND_MAX * L) - Lmoitie;
        p.y = ((double)rand() / RAND_MAX * L) - Lmoitie;

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
int change_nouvelle_configuration(Particule tab_par[], int nbx_particules, double L_moitie, double Lboite)
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
    if (tab_par[index_particule].x >= L_moitie) {tab_par[index_particule].x -=Lboite;}
    if (tab_par[index_particule].x < -L_moitie) {tab_par[index_particule].x += Lboite;}

    //Application des conditions périodiques sur y : 
    if (tab_par[index_particule].y >= L_moitie) {tab_par[index_particule].y -=Lboite;}
    if (tab_par[index_particule].y < -L_moitie) {tab_par[index_particule].y += Lboite;}

    return index_particule;
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
double calcul_energie_potentielle_terme_commun(Particule tab_par[], int index_modifie) 
{
    /*définition des variables locales à utiliser dans les calculs*/
    double U = 0.0; // energie potentielle
    double rij2, inv_r2, inv_r2x3, inv_r2x6; //puissance de rij
    double V = L * L; //volume

    for (int i = 0; i < nbx_particules; i++) 
    {

        if (i!=index_modifie){

            /*condition qui verifie que le calcul se fait que par rapports aux partiucules et ignore les défauts*/
            if (tab_par[i].actif) 
            {
                for (int j = i + 1; j < nbx_particules; j++) 
                {

                    if (j!=index_modifie){

                        if (!tab_par[j].actif) continue; 

                        double x_sous = tab_par[j].x - tab_par[i].x; //Dx
                        double y_sous = tab_par[j].y - tab_par[i].y; //

                        double abs_distx = fabs(x_sous);
                        double abs_disty = fabs(y_sous);

                        if (abs_distx > Lmoitie) {
                            x_sous = - (x_sous / abs_distx) * (L - abs_distx);
                        }

                        if (abs_disty > Lmoitie) {
                            y_sous = - (y_sous / abs_disty) * (L - abs_disty);
                        }

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
        }
        
    }
    return U;    //calcul et mise à jour l'energie potentielle sans U_tail
}


double calcul_energie_potentielle_terme_modifie(Particule tab_par[], int index_modifie) 
{
    /*définition des variables locales à utiliser dans les calculs*/
    double U_index = 0.0; // energie potentielle
    double rij2, inv_r2, inv_r2x3, inv_r2x6; //puissance de rij
    double V = L* L; //volume
    int dimension = 2;

    int i = index_modifie; 

    if (tab_par[i].actif){

        for (int j=0; j<nbx_particules; j++){
            if (!tab_par[j].actif) continue; 

            if(j!=i){
                double x_sous = tab_par[j].x - tab_par[i].x; //Dx
                double y_sous = tab_par[j].y - tab_par[i].y; //

                double abs_distx = fabs(x_sous);
                double abs_disty = fabs(y_sous);

                if (abs_distx > Lmoitie) {
                    x_sous = - (x_sous / abs_distx) * (L - abs_distx);
                }

                if (abs_disty > Lmoitie) {
                    y_sous = - (y_sous / abs_disty) * (L - abs_disty);
                }

                /*Calcul de la distance entre deux particules*/
                rij2 = x_sous * x_sous + y_sous * y_sous;   

                /*Calcul sur les puissances de rij pour eviter l'appel de pow*/
                inv_r2 = 1.0 / rij2;                  //   1/rij^2
                inv_r2x3 = inv_r2 * inv_r2 * inv_r2; //    1/rji^6
                inv_r2x6 = inv_r2x3 * inv_r2x3;     //    1/rij^12

                /*Conditions qui permet de tronquer le potentiel à Rc*/
                if (rij2 < Rc2) 
                {
                    U_index += 4 * (inv_r2x6 - inv_r2x3) - U_decal; //calcul du potentiel entre deux particules
                }
            }
        }
    }

    return U_index;
}

//-----------------------------{calcul de boltzmann}------------------------------------------

double calcul_boltzmann(Particule tab_ini[], Particule tab_temp[], int index_modifie)
{
    double proba = 0.0;

    double terme_commun = calcul_energie_potentielle_terme_commun(tab_ini, index_modifie) + U_tail;

    double energie_tab_ini = terme_commun + calcul_energie_potentielle_terme_modifie(tab_ini,index_modifie); 
    double energie_tab_temp = terme_commun + calcul_energie_potentielle_terme_modifie(tab_temp,index_modifie); 

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


//-----------------------------{main}----------------------------------------------------------------
int main() 
{

    // Initialiser le générateur de nombres aléatoires avec une graine pour avoir des résultats reproductibles
    srand(time(NULL));

    // Creation d'un tableau de particules
    Particule tab_particules[nbx_particules];
    Particule tab_particules_temporaire[nbx_particules]; 


    // Initialisation des particules avec une configuration aleatoire
    initialiser_Configuration_aleatoire(tab_particules);

    printf("aaa \n");
    
    // Ouvre le fichier "positions_data.txt" en mode écriture ("w")
    // Ce fichier stockera les positions des particules au fil du temps dans le format suivant:
    // temps   x1 y1    x2 y2    x3 y3    x4 y4.....
    FILE *pos_file = fopen("positions_data.txt", "w");

    int tentative = 0;
    int tentative_acceptee=0; 
    int index_modidi=0; 

    // Appelle la fonction enregistrer_Positions pour écrire les positions initiales des particules
    enregistrer_positions(pos_file, tab_particules, tentative, nbx_particules);

    for (int cycle=0 ; cycle<nombre_cycles; cycle++){

        printf("Cycle : %d\n", cycle); 

            for (int tentativou=0; tentativou < nbx_particules ; tentativou++){

                for (int i=0; i<nbx_particules; i++){
                    tab_particules_temporaire[i]=tab_particules[i];
                }

                int index_modidi = change_nouvelle_configuration(tab_particules_temporaire, nbx_particules, Lmoitie, L); 
                double proba_finale=calcul_boltzmann(tab_particules, tab_particules_temporaire,index_modidi); 

                double x= ((double)rand() / RAND_MAX);

                if (proba_finale >= x){

                    tentative_acceptee += 1;

                    for (int i=0; i<nbx_particules; i++){
                        tab_particules[i]=tab_particules_temporaire[i];
                    }
                }

                tentative += 1; 

            }

        enregistrer_positions(pos_file, tab_particules, tentative, nbx_particules); 
    }

    fclose(pos_file);  // De même pour le fichier positions_data.txt

    double taux_acceptation = (double)tentative_acceptee/tentative; 
    printf("Taux d'acceptation : %f", taux_acceptation); 

    return 0;
}
