//////////
//HEADER//
//////////
#ifndef HEADER_H
#define HEADER_H

#define epsilon_0 = 8.85e-1

////////////
//INCLUDES//
////////////

#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <vector>
#include <functional>
#include <cstring>

//////////////
//STRUCTURES//
//////////////

struct position
{
  double x;
  double y;
  double z;
};

struct vecteur
{
  position ptApp;

  double compoX
  double compoY;
  double compoZ;

  double norme;
};

struct force
{
  vecteur vecForce;
  double valeur;
  position lieuApp;
}

struct partChargee
{
  double masse;
  double charge;

  position lieu;
  vecteur vitesse;
};

////////////////////////////////////////////////////
//FONCTIONS POUR L'INITIALISATION DE LA SIMULATION//
////////////////////////////////////////////////////

void recuperation_donnees(double&, int&, std::string&, double&, int&, const std::string &);
int methode_recup_donnees();
int choix_config();
std::vector <std::string> ReadInputsFromFile(const std::string& filename);
void recuperation_donnees_fichier(double &, int &, std::string &, double &, int &, const std::string&);
void recuperation_donnees_console(double &, int &, std::string &, double &, int &);

////////////////////////
//PLACEMENT DES ATOMES//
////////////////////////

//Fonctions intermédiaires (formes)

void placer_cercle(int &, double &, double* &);
void placer_carre(int &, double &, double* &);
void placer_triangle(int &, double &, double* &);
void placement_aleatoire(int&, double&, double*&);
void placer_sphere(int&, double&, double*&);


//Fonction principale

void placer_atomes(std::string &, int &, double &, double* &);
void placer_atomes_3D(std::string &, int &, double &, double* &);

//Calcul de la distance

double calcul_distance_3D(double, double, double);
double calculDistance(position, position);
double calcul_potentiel_3D(int, int, double*&);


///////////////////////
//FONCTIONS DE CALCUL//
///////////////////////

double calcul_NRJ_cov(double*&, const int &);
double calcul_NRJ_vdw(double*&, const int &);
double calcul_piege(int&, double*&);
double calcul_vdw(int&, int&, double*&);
double calcul_cov(int&, int&, double*&);

void Minimiser_Vdw(double* &, const int&, double&, const double&, const int &);
void Minimiser_Cov(double* &, const int&, double&, const double&, const int &);

force calcForceElec();

void affichage_vecteur(double*&, const int &);
void affichage_vecteur_3D(double*&, const int &);

//SURCHARGE OPERATEURS ETC//
vecteur operator-(const vecteur&, const vecteur&);
vecteur operator+(const vecteur&, const vecteur&);
vecteur operator*(const vecteur&, const vecteur&);
vecteur operator*(const double&, const vecteur&);

#endif