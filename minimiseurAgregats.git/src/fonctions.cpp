/////////////
//Fonctions//
/////////////

/////////
//LINKS//
/////////
#include "header.h"

////////////
//INCLUDES//
////////////
#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <cstring>

//using namespace std ;

int methode_recup_donnees()
{
  int numero_methode = 0;
  
  std::cout << "Options pour initialiser la configuration de depart : " << std::endl;
  std::cout << "1. A partir d'un fichier texte (recommande)" << std::endl;
  std::cout<< "2. Par interaction avec la console" << std::endl;
  std::cout<< "3. Configuration par defaut" << std::endl << std::endl;
  std::cout<< "Entrez le numero de la methode choisie svp" << std::endl;

  std::cin >> numero_methode ;

  return numero_methode;
}

void recuperation_donnees(double& mu, int& n_atomes, std::string& type_forme, double& param_supp, int& type_potentiel, const std::string & filename)
{
  int n = methode_recup_donnees();

  switch(n)
    {
      case 1 :
      recuperation_donnees_fichier(mu, n_atomes, type_forme, param_supp, type_potentiel, filename) ;
      std::cout << "Donnees recuperees a partir du fichier " << filename << std::endl;
      break ;

      case 2 : 
      recuperation_donnees_console(mu, n_atomes, type_forme, param_supp, type_potentiel);
      std::cout << "Donnees recuperees " << filename << std::endl;
      break;

      case 3 :
      std::cout << "La configuration par defaut sera utilisee" << std::endl;
      break;

      default :
      std::cerr << "Le numero entre n'est pas valide" << std::endl ;
    }
  
  
}

int choix_config()
{
  int numero_config ;
  std::cout << "Numero configuration souhaitee ?" << std::endl ;
  std::cin >> numero_config ;
  return numero_config ;
}

std::vector <std::string> ReadInputsFromFile(const std::string& filename)
{
  std::vector <std::string> results;
  std::ifstream file(filename);
  std::string buffer;
  if(file.is_open())
  {
    while(getline(file, buffer))
      {
        results.push_back(buffer);
      }
    file.close();
  }

  else
  {
    std::cerr << "Unable to open file" << filename << std::endl ;
  }
  return results;
}

void recuperation_donnees_fichier(double & mu, int & n_atomes, std::string & type_forme, double & param_supp, int & type_potentiel, const std::string& filename)
{
  int i = choix_config() ;

  mu = stod(ReadInputsFromFile(filename)[(6*i)+1], 0);
  n_atomes = stoi(ReadInputsFromFile(filename)[(6*i)+2], 0);
  type_forme = ReadInputsFromFile(filename)[(6*i)+3];
  param_supp = stod(ReadInputsFromFile(filename)[(6*i)+4], 0);
  type_potentiel = 0;
}

void recuperation_donnees_console(double& mu, int& n_atomes, std::string& type_forme, double& param_supp, int& type_potentiel)
{
  int tf = 0 ;
  
  std::cout << "Profondeur du piege souhaitee ? " << std::endl;
  std::cin >> mu;
  
  std::cout << "Nombre d'atomes souhaite? " << std::endl;
  std::cin >> n_atomes;
  
  std::cout << "Forme pour la disposition des atomes ? (choisir un numero parmi les options proposees) " << std::endl;
  std::cout << "1. Cercle" << std::endl;
  std::cout << "2. Carre" << std::endl;
  std::cout << "3. Triangle" << std::endl;
  std::cout << "4. Aleatoire" << std::endl;
  
  std::cin >> tf;
    switch(tf)
      {
        case 1 :
          type_forme = "Cercle" ;
          break;

        case 2 : 
          type_forme = "Carre" ;
          break ;

        case 3 : 
          type_forme = "Triangle" ;
          break ;

        case 4 :
          type_forme = "Aleatoire" ;
          break ;
        }
  
  std::cout << "Type de potentiel ? (choisir un numero parmi les options proposees)" << std::endl;
  std::cout << "1. Van der Waals" << std::endl;
  std::cout << "2. Covalent" << std::endl ;
  std::cin >> type_potentiel;
}

////////////////////////
//Placement des atomes//
////////////////////////

//Fontions intermédiaires

void placer_cercle(int & n_atomes, double & param_supp, double* & P_coordonnees)
{
  for(int i=0; i<n_atomes ; i++)
    {
      double x = param_supp*cos((2.0*(i+1)*M_PI)/n_atomes) ;
      double y = param_supp*sin((2.0*(i+1)*M_PI)/n_atomes);
      P_coordonnees[2*i] = x ;
      P_coordonnees[(2*i)+1] = y ;
    }
}

void placer_carre(int & n_atomes, double & param_supp, double* & P_coordonnees)
{
  for(int i=0; i<n_atomes ; i++)
    {
      double x_max = param_supp/2;
      double x_min = -(param_supp/2);
      double y_max = param_supp/2;
      double y_min = -(param_supp/2);


      //Remplissage sommet 1 à 2
      for(int j=0 ; j<n_atomes/4; j++)
        {
          double x = x_min + (j*param_supp/(n_atomes/4)) ;
          double y = y_max ;
          P_coordonnees[2*i] = x ;
          P_coordonnees[(2*i)+1] = y ;
          i++;
        }

      //Remplissage sommet 2 à 3
      for(int j=0; j<n_atomes/4; j++)
        {
          double x = x_max ;
          double y = y_max - (j*param_supp/(n_atomes/4)) ;
          P_coordonnees[2*i] = x ;
          P_coordonnees[(2*i)+1] = y ;
          i++;
        }

      //Remplissage sommet 3 à 4
      for(int j=0; j<n_atomes/4; j++)
        {
          double x = x_max - (j*param_supp/(n_atomes/4)) ;
          double y = y_min ;
          P_coordonnees[2*i] = x ;
          P_coordonnees[(2*i)+1] = y ;
          i++;
        }

      //Remplissage sommet 4 à 1
      for(int j=0; j<n_atomes/4; j++)
        {
          double x = x_min ;
          double y = y_min + (j*param_supp/(n_atomes/4));
          P_coordonnees[2*i] = x ;
          P_coordonnees[(2*i)+1] = y ;
          i++;
        }
    }
}

void placer_triangle(int & n_atomes, double & param_supp, double* & P_coordonnees)
{
  double hauteur_triangle = (sqrt(3.0)/2.0)*param_supp ;
  int n_atomes_par_cote = n_atomes/3 ;
  int j = 0 ;

  //Sommets : 
  double sommet_1_x = 0.0 ;
  double sommet_1_y = (2.0/3.0)*hauteur_triangle ;

  double sommet_2_x = param_supp/2.0 ;
  double sommet_2_y = -(1.0/3.0)*hauteur_triangle ;

  double sommet_3_x = -(param_supp/2.0) ;
  double sommet_3_y = -(1.0/3.0)*hauteur_triangle ;

  //Remplissage sommet 1 à 2

  for(int i=0; i<n_atomes_par_cote ; i++)
    {
      double k = 0.0 ;
      double x = 0.0 ;
      double y = 0.0 ;

      k = static_cast <double> (i)/(n_atomes_par_cote - 1.0) ;
      x = ((1.0-k) * sommet_1_x) + (k * sommet_2_x) ;
      y = ((1.0-k) * sommet_1_y) + (k * sommet_2_y) ;

      P_coordonnees[2*j] = x ;
      P_coordonnees[1 + 2*j] = y ;
      j++ ;
    }

  //Remplissage sommet 2 à 3

  for(int i=0; i<n_atomes_par_cote ; i++)
    {
      double k = 0.0 ;
      double x = 0.0 ;
      double y = 0.0 ;

      k = static_cast <double> (i)/(n_atomes_par_cote - 1.0) ;
      x = ((1.0-k) * sommet_2_x) + (k * sommet_3_x) ;
      y = ((1.0-k) * sommet_2_y) + (k * sommet_3_y) ;

      P_coordonnees[2*j] = x ;
      P_coordonnees[1 + 2*j] = y ;
      j++ ;
    }

  //Remplissage sommet 3 à 1

  for(int i=0; i<n_atomes_par_cote ; i++)
    {
      double k = 0.0 ;
      double x = 0.0 ;
      double y = 0.0 ;

      k = static_cast <double> (i)/(n_atomes_par_cote - 1.0) ;
      x = ((1.0-k) * sommet_3_x) + (k * sommet_1_x) ;
      y = ((1.0-k) * sommet_3_y) + (k * sommet_1_y) ;

      P_coordonnees[2*j] = x ;
      P_coordonnees[1+ 2*j] = y ;
      j++ ;
    }
}

void placement_aleatoire(int & n_atomes, double & param_supp, double* & P_coordonnees)
{
  
  srand(time(0));
  for (int i = 0 ; i < n_atomes ; i ++)
    {

      double x = 0.0 ;
      double y = 0.0 ;


      x = ((1.0*rand())/RAND_MAX)*param_supp*((2.0*(1.0*rand()/RAND_MAX)) - 1.0) ;

      y = ((1.0*rand())/RAND_MAX)*param_supp*((2.0*(1.0*rand()/RAND_MAX)) - 1.0) ;

      P_coordonnees[2*i] = x ;
      P_coordonnees[(2*i)+1] = y ;
    }
}

void placer_sphere(int& n_atomes, double& param_supp, double*& P_coordonnees/*, int& maxiter_placement, double& pas_placement*/)
{
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
  double theta = 0.0;
  double rayonBis = 0.0;

  double goldenRatio = (1.0 + sqrt(5.0))/2.0 ;

//a) Fibonacci (a venir)

  for(int i = 0; i<n_atomes; i++)
    {
      rayonBis = sqrt(1.0 - z * z);
      theta = (i*2.0*M_PI)/goldenRatio ;
      x = param_supp*rayonBis*cos(theta);
      y = param_supp*rayonBis*sin(theta);
      z = 1.0 -((2.0*i)/(n_atomes-1.0));

      P_coordonnees[3*i] = x;
      P_coordonnees[3*i+1] = y;
      P_coordonnees[3*i+2] = z;
    }
      

//b) Minimisation du modele de Thomson

  //On effectue un premier placement simple, dans un plan
/*  srand(time(0));
  double U_P_coordonnees = 0.0 ;
  double U_nouveau_vec = 0.0 ;
  double valeur_deplacement = 0.0 ;
  double* nouveau_vec = new double[3*n_atomes] ;
  
  for(int i=0; i<n_atomes ; i++)
    {
      double x = param_supp*cos((2.0*(i+1)*M_PI)/n_atomes) ;
      double y = param_supp*sin((2.0*(i+1)*M_PI)/n_atomes);
      P_coordonnees[3*i] = x ;
      P_coordonnees[(3*i)+1] = y ;
      P_coordonnees[(3*i)+2] = 0 ; //Correspond à la coordonnée z
    }*/


//On fait ensuite une minimisation
/*  for(int i = 0 ; i < maxiter_placement ; i ++)
    {
      for (int j = 0 ; j<n_atomes*3 ; j ++)
        {
          nouveau_vec[j] = P_coordonnees[j] ;
        }

      for(int j = 0 ; j<n_atomes*3 ; j ++)
        {
          valeur_deplacement = (2.0*(1.0*rand()/RAND_MAX)) - 1.0 ;
          valeur_deplacement *= pas_placement ;
          nouveau_vec[j] += valeur_deplacement ;
        }

      for(int k = 0; k<n_atomes*3 ; k++)
        {
          for (int l = 1; l < n_atomes*3 ; l++)
            {
              U_nouveau_vec += calcul_potentiel_3D(k, l, nouveau_vec);
              U_P_coordonnees += calcul_potentiel_3D(k, l, P_coordonnees);
            }

      if(U_nouveau_vec < U_P_coordonnees)
      {
        for(int j = 0 ; j < n_atomes*3 ; j++)
          {
            P_coordonnees[j] = nouveau_vec[j] ;
          }
      }
    }
  
  }*/
}

/*double calcul_potentiel_3D(int a, int b, double*& vec)
{
  double result = 0.0 ;
  double x_i = 0.0 ;
  double y_i = 0.0 ;
  double z_i = 0.0 ;
  double x_j = 0.0 ;
  double y_j = 0.0 ;
  double z_j = 0.0 ;

  x_i = vec[3*a] ;
  y_i = vec[1+3*a] ;
  z_i = vec[2+3*a] ;
  x_j = vec[3*b] ;
  y_j = vec[1+3*b] ;
  z_j = vec[2+3*b] ;

  double d_x = x_i - x_j ;
  double d_y = y_i - y_j ;
  double d_z = z_i - z_j ;

  double distance = calcul_distance_3D(d_x, d_y, d_z) ;

  result = 1.0/distance ;
  return result ;
}*/


//Fonction principale

void placer_atomes(std::string & type_forme, int & n_atomes, double & param_supp, double* & P_coordonnees)
{
  if(type_forme == "Cercle")
  {
    placer_cercle(n_atomes, param_supp, P_coordonnees);
  }

  else if(type_forme == "Carre")
  {
    placer_carre(n_atomes, param_supp, P_coordonnees);
  }

  else if(type_forme == "Triangle_Equilateral")
  {
    placer_triangle(n_atomes, param_supp, P_coordonnees);
  }

  else if(type_forme == "Aleatoire")
  {
    placement_aleatoire(n_atomes, param_supp, P_coordonnees);
  }

  else 
  {
    std::cerr << "Erreur : forme non reconnue" << std::endl ; ;
  }

}

void placer_atomes_3D(std::string & type_forme, int & n_atomes, double & param_supp, double* & P_coordonnes)
{
  return ;
}


/////////////////////
//TEST PLACEMENT 3D//
/////////////////////

/*
struct structPlacement
{
  double q;
  double x;
  double y;
  double z;
  double m;
};

void test3D(int n)
{
  std::vector
  for(int i = 0; i<n; i++)
  (

    Minimisation

  )

}

vecteur deplacement

calcul force elec

*/


void affichage_vecteur(double*& P_coordonnees, const int & n_atomes)
{
  for(int i = 0 ; i < n_atomes ; i++)
    {
      std::cout << std::setw(15) << P_coordonnees[i*2] << " " << P_coordonnees[i*2 + 1] << " ";
    }
  std::cout << std::endl ;
}

void affichage_vecteur_3D(double*& P_coordonnees, const int & n_atomes)
{
  for(int i = 0 ; i < n_atomes ; i++)
    {
      std::cout << std::setw(15) << P_coordonnees[i*3] << " " << P_coordonnees[i*3 + 1] << " " << P_coordonnees[i*3 + 2] << " ";
    }
  std::cout << std::endl ;
}
