///////////////////////////////////////////
//MAIN - MINIMISEUR D'UN AGREGAT D'ATOMES//
///////////////////////////////////////////

/////////
//LINKS//
/////////

#include "Header.h"

////////////
//INCLUDES//
////////////

#include <cstring> //memcpy
#include <fstream>
#include <iostream>
#include <iomanip>

//using namespace std ;

////////////////////////////
//DECLARATION DE VARIABLES//
////////////////////////////

int n_atomes = 20 ;
int double_n_atomes = 40 ;
int triple_n_atomes = 60 ;
double profondeur_mu = 0.0 ; //"Profondeur" mu du piège ; plus elle est élevée, plus l'énergie potentielle augmente
//Lorsque l'on s'éloigne du centre

std::string type_forme = "Cercle "; //Stockage de la forme de départ dans laquelle sont disposés les atomes à l'instant t=0 de 
//la simulation
double param_supp = 0.0 ;
int type_potentiel = 1; //Type de potentiel que l'on souhaite calculer (1 pour Van der Waals, 2 pour covalent)
std::string filename = "Donnees.txt" ; //Fichier à partir duquel on peut récupérer des données de départ. Cf "..." pour
//la configuration du fichier

double pas = 0.0 ;
double decremente = 0.0 ;
int maxiter = 0 ; //Nombre max. d'itérations = combien de fois on essaye de minimiser l'énergie avant la configuration
//finale

//Mesure des temps mis pour les différents potentiels
double temps_vdw = 0.0 ;
double dt_vdw = 0.0 ;
double temps_cov = 0.0 ;
double dt_cov = 0.0 ;

//3D feature
int numberOfDimensions = 0;


/////////////////
//CORPS DU MAIN//
/////////////////

int main( int argc, char **argv )
{
    //Recuperation de toutes les donnees contenues dans le fichier "Donnees.txt"
    ReadInputsFromFile(filename);

    //Choix d'une configuration particulière du fichier "Donnees.txt" par l'utilisateur
    recuperation_donnees(profondeur_mu, n_atomes, type_forme, param_supp, type_potentiel, filename);
    double_n_atomes = n_atomes*2 ;
    triple_n_atomes = n_atomes*3 ;

    //Affichage des donnees dans le terminal
    std::cout << "Parametres configuration : " << std::endl ;
    std::cout << "Piege avec une profondeur mu = " << profondeur_mu << std::endl ;
    std::cout << "Contient " << n_atomes << " atomes disposes en " << type_forme << "." << std::endl ;

    std::cout << "Valeur du pas souhaitee ?" << std::endl ;
    std::cin >> pas ;
    std::cout << "Valeur du facteur de decrementation souhaitée ?" << std::endl ; 
    std::cin >> decremente ;
    std::cout << "Valeur de l'iteration maximum souhaitee ?" << std::endl ;
    std::cin >> maxiter ;
    std::cout << "Type de potentiel souhaite ? " << std::endl ;
    std::cout << "1. Van der Waals" << std::endl ;
    std::cout << "2. Covalent" << std::endl ;
    std::cin >> type_potentiel ;

    double* P_coordonnees_1 = new double[double_n_atomes];
    double* P_coordonnees_2 = new double[double_n_atomes];
    double* P_coordonnees_1_sphere = new double[triple_n_atomes];
    double* P_coordonnees_2_sphere = new double[triple_n_atomes];

    placer_atomes(type_forme, n_atomes, param_supp, P_coordonnees_1) ;
    placer_atomes(type_forme, n_atomes, param_supp, P_coordonnees_2) ;

    if(type_potentiel == 1)
        {
            temps_vdw = clock ();    
            Minimiser_Vdw(P_coordonnees_1, n_atomes, pas, decremente, maxiter) ;
            dt_vdw = clock() - temps_vdw ;
            dt_vdw /= (double)CLOCKS_PER_SEC ;
        
            std::cout << "Minimisation pour Van der Waals terminée, avec un temps de calcul de " << dt_vdw << " secondes." << std::endl ;
        }

    else if(type_potentiel == 2)
        {
            temps_cov = clock ();    
            Minimiser_Cov(P_coordonnees_2, n_atomes, pas, decremente, maxiter) ;
            dt_cov = clock() - temps_cov ;
            dt_cov /= (double)CLOCKS_PER_SEC ;
        
            std::cout << "Minimisation pour Covalent terminée, avec un temps de calcul de " << dt_cov << " secondes." << std::endl ;
        }

    else
        {
            temps_vdw = clock ();    
            Minimiser_Vdw(P_coordonnees_1, n_atomes, pas, decremente, maxiter) ;
            dt_vdw = clock() - temps_vdw ;
            dt_vdw /= (double)CLOCKS_PER_SEC ;
        
            std::cout << "Minimisation pour Van der Waals terminée, avec un temps de calcul de " << dt_vdw << " secondes." << std::endl ;

            temps_cov = clock ();    
            Minimiser_Cov(P_coordonnees_2, n_atomes, pas, decremente, maxiter) ;
            dt_cov = clock() - temps_cov ;
            dt_cov /= (double)CLOCKS_PER_SEC ;
        
            std::cout << "Minimisation pour Covalent terminée, avec un temps de calcul de " << dt_cov << " secondes." << std::endl ;
        }

    return 0 ; 
}


