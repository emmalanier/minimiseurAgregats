/////////////
//Fonctions//
/////////////

/////////
//LINKS//
/////////
#include "Header.h"

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

void placer_sphere(int& n_atomes, double& param_supp, double*& P_coordonnees)
{
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;

//a) Fibonacci

//b) Minimisation du modele de Thomson

  //On effectue un premier placement simple, dans un plan
    for(int i=0; i<n_atomes ; i++)
    {
      double x = param_supp*cos((2.0*(i+1)*M_PI)/n_atomes) ;
      double y = param_supp*sin((2.0*(i+1)*M_PI)/n_atomes);
      P_coordonnees[3*i] = x ;
      P_coordonnees[(3*i)+1] = y ;
      P_coordonnees[(3*i)+2] = 0 ; //Correspond à la coordonnée z
    }

  //On fait ensuite une minimisation
  
}

/*bool plus_ou_moins()
{
//Fonction temporaire
  
  bool pm = false ;
  double pm_val = 0.0 ;

  srand(time(0)) ;
  pm_val = (1.0*rand()/RAND_MAX) ;

  if (pm_val < 0.5)
    pm = true ;

  else 
    pm = false ;

  return pm ;
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

}

double calculer_distance_3D(double x, double y, double z)
{
  double results = 0.0;
  double results_inter = 0.0;

  results_inter = (x*x)+(y*y)+(z*z);

  results = 
}

//Energie et potentiel//

double calcul_cov(int& a, int& b, double*& vec) //Calcule le potentiel covalent entre un atome i et un atome j
{
  double result = 0.0 ;
  double x_i = 0.0 ; //Var pour stock coordonnee en x de l'atome i
  double y_i = 0.0 ; //Var pour stock coordonnee en y de l'atome i
  double x_j = 0.0 ; //Var pour stock coordonnee en x de l'atome j
  double y_j = 0.0 ; //Var pour stock coordonnee en y de l'atome j

  //Transmission des coordonnees aux variables prévues pour
  x_i = vec[2*a] ;
  y_i = vec[1+2*a] ;
  x_j = vec[2*b] ;
  y_j = vec[1+2*b] ;

  //Calcul des distances entre les coordonnees des deux atomes
  double d_x = x_i - x_j ;
  double d_y = y_i - y_j ;

  double r_carre = (d_x * d_x) + (d_y * d_y) ; //Distance entre les deux atomes au carré
  double u_carre = r_carre + 0.25 ;
  double s_carre = (d_x * d_x) - (d_y * d_y) ;
  double s_puissance_4 = s_carre * s_carre ;
  double u_puissance_4 = u_carre * u_carre ;

  result =  (1.0/u_puissance_4) * (1.0 - (s_puissance_4/u_carre)); //Formule donnant la valeur du potentiel
  return result ;
}

double calcul_vdw(int& a, int& b, double*& vec)
{
  double result = 0.0 ;
  double x_i = 0.0 ;
  double y_i = 0.0 ;
  double x_j = 0.0 ;
  double y_j = 0.0 ;

  x_i = vec[2*a] ;
  y_i = vec[1+2*a] ;
  x_j = vec[2*b] ;
  y_j = vec[1+2*b] ;

  double d_x = x_i - x_j ;
  double d_y = y_i - y_j ;

  double r_carre = (d_x * d_x) + (d_y * d_y) ;
  double u_carre = r_carre + 0.25 ;
  double u_puissance_4 = u_carre * u_carre ;

  result = (1.0/u_puissance_4) - (1.0/u_carre) ;
  return result ;
}

double calcul_piege(int& a, double*& vec) //Calcule l'énergie potentielle du piège en 1 point
{
  double coordo_x = 0.0 ;
  double coordo_y = 0.0 ;
  double result = 0.0 ;

  coordo_x = vec[2*a] ;
  coordo_y = vec[1+2*a] ;

  double x_carre = coordo_x * coordo_x ;
  double y_carre = coordo_y * coordo_y ;

  result = 0.1 * (x_carre + y_carre) ;

  return result ;
}

double calcul_NRJ_vdw(double*& vec, const int & n_atomes) //Somme les potentiels entre chaque atomes, pour obtenir
//L'énergie totale du système
{
  double result = 0.0;

  for(int i = 0; i < n_atomes; i++) 
    {
      result += calcul_piege(i, vec);
      for(int j = i + 1; j < n_atomes; j++) 
        {
          result += calcul_vdw(i, j, vec);
        }
    }

  return result;
}

double calcul_NRJ_cov(double*& vec, const int & n_atomes)
{
  double result = 0.0;

  for(int i = 0; i < n_atomes; i++) 
    {
      result += calcul_piege(i, vec);
      for(int j = i + 1; j < n_atomes; j++) 
        {
          result += calcul_cov(i, j, vec);
        }
    }

  return result;
}


/////////////////////////////
//FONCTIONS DE MINIMISATION//
/////////////////////////////

//Pour minimiser les potentiels, ici, on passe par une recherche de Monte-Carlo. L'algorithme, dans la première 
//partie, va créer un nouveau vecteur, identique à celui stockant les coordonnées des atomes. Dans ce nouveau vecteur, 
//les coordonnées vont TOUTES être modifiées de manière aléatoire. L'énergie correspondant au système que représente 
//le nouveau vecteur est alors calculée. Si celle-ci est inférieure à celle du vecteur comportant les coordonnées 
//actuelles des atomes, les valeurs du nouveau vecteur remplacent les anciennes. Sinon, le vecteur des coordonnées
//actuelles reste inchangé.

//Une fois que l'on a fait un quart des itérations, le fonctionnement diffère légèrement. Au lieu de changer les 
//coordonnées de tous les atomes puis d'évaluer l'énergie du nouveau système, seules les coordonnées d'un atome sont 
//modifiées avant de réévaluer l'énergie. Si cette modification permet de la diminuer, elle est conservée. Sinon, 
//l'atome conserve ses coordonnées actuelles.

void Minimiser_Vdw(double* & P_coordo, const int& n_atomes, double& pas, const double & decremente, const int & maxiter)
{
  std::streambuf* coutbuf = std::cout.rdbuf();
  std::ofstream fichier_log_NRJ("NRJ_minimisation_vdw.txt");
  std::ofstream fichier_log_positions("Donnees_minimisation_vdw.txt");

  srand(time(0));
  double NRJ_P_coordo = 0.0 ;
  double NRJ_nouveau_vec = 0.0 ;
  double valeur_deplacement = 0.0 ;
  double* nouveau_vec = new double[2*n_atomes] ;


  for(int i = 0 ; i < maxiter ; i ++)
    {
      //Partie 1 de la recherche de Monte-Carlo : tout le vecteur est modifié avant d'évaluer l'énergie
      for (int j = 0 ; j<n_atomes*2 ; j ++)
        {
          nouveau_vec[j] = P_coordo[j] ;
        }

      for(int j = 0 ; j<n_atomes*2 ; j ++)
        {
          valeur_deplacement = (2.0*(1.0*rand()/RAND_MAX)) - 1.0 ;
          valeur_deplacement *= pas ;
          nouveau_vec[j] += valeur_deplacement ;
        }

      NRJ_P_coordo = calcul_NRJ_vdw(P_coordo, n_atomes);
      NRJ_nouveau_vec = calcul_NRJ_vdw(nouveau_vec, n_atomes);

      if(NRJ_nouveau_vec < NRJ_P_coordo)
      {
        for(int j = 0 ; j < n_atomes*2 ; j++)
          {
            P_coordo[j] = nouveau_vec[j] ;
          }
      }

      //Partie 2 de la recherche de Monte-Carlo : 1 modification de coordonnées = 1 évaluation de l'énergie
      if(i >= maxiter/4)
      {
        for(int j = 0 ; j < n_atomes ; j++)
          {
            valeur_deplacement = (2.0*(1.0*rand()/RAND_MAX)) - 1.0 ;
            valeur_deplacement *= pas ;
            nouveau_vec[j] += valeur_deplacement ;

            NRJ_P_coordo = calcul_NRJ_vdw(P_coordo, n_atomes);
            NRJ_nouveau_vec = calcul_NRJ_vdw(nouveau_vec, n_atomes);  

            if(NRJ_nouveau_vec < NRJ_P_coordo)
            {
              P_coordo[j] = nouveau_vec[j];
            }
          }
      }

      if(i%20 == 0) //A modifier pour que l'utilisateur puisse moduler la fréquence d'output des données
      {
        std::cout.rdbuf(fichier_log_NRJ.rdbuf());
        std::cout << "Iter #" << i << " :" << std::endl ;
        std::cout << "Ancienne energie : " << NRJ_P_coordo << std::endl ;
        std::cout << "Nouvelle energie : " << NRJ_nouveau_vec << std::endl ;
        std::cout.rdbuf(fichier_log_positions.rdbuf());
        affichage_vecteur(P_coordo, n_atomes) ;
        pas/=decremente;
      }

    }
  std::cout.rdbuf(coutbuf);
  fichier_log_NRJ.close();
  fichier_log_positions.close();  
}

void Minimiser_Cov(double* & P_coordo, const int& n_atomes, double& pas, const double & decremente, const int & maxiter)
{
  std::streambuf* coutbuf = std::cout.rdbuf();
  std::ofstream fichier_log_NRJ("NRJ_minimisation_cov.txt");
  std::ofstream fichier_log_positions("Donnees_minimisation_cov.txt");

  srand(time(0));
  double NRJ_P_coordo = 0.0 ;
  double NRJ_nouveau_vec = 0.0 ;
  double valeur_deplacement = 0.0 ;
  double* nouveau_vec = new double[2*n_atomes] ;


  for(int i = 0 ; i < maxiter ; i ++)
    {
      for (int j = 0 ; j<n_atomes*2 ; j ++)
        {
          nouveau_vec[j] = P_coordo[j] ;
        }

      for(int j = 0 ; j<n_atomes*2 ; j ++)
        {
          valeur_deplacement = (2.0*(1.0*rand()/RAND_MAX)) - 1.0 ;
          valeur_deplacement *= pas ;
          nouveau_vec[j] += valeur_deplacement ;
        }

      NRJ_P_coordo = calcul_NRJ_cov(P_coordo, n_atomes);
      NRJ_nouveau_vec = calcul_NRJ_cov(nouveau_vec, n_atomes);

      if(NRJ_nouveau_vec < NRJ_P_coordo)
      {
        for(int j = 0 ; j < n_atomes*2 ; j++)
          {
            P_coordo[j] = nouveau_vec[j] ;
          }
      }

      if(i >= maxiter/4)
      {
        for(int j = 0 ; j < n_atomes ; j++)
          {
            valeur_deplacement = (2.0*(1.0*rand()/RAND_MAX)) - 1.0 ;
            valeur_deplacement *= pas ;
            nouveau_vec[j] += valeur_deplacement ;

            NRJ_P_coordo = calcul_NRJ_cov(P_coordo, n_atomes);
            NRJ_nouveau_vec = calcul_NRJ_cov(nouveau_vec, n_atomes);  

            if(NRJ_nouveau_vec < NRJ_P_coordo)
            {
              P_coordo[j] = nouveau_vec[j];
            }
          }
      }

      if(i%20 == 0)
      {
        std::cout.rdbuf(fichier_log_NRJ.rdbuf());
        std::cout << "Iter #" << i << " :" << std::endl ;
        std::cout << "Ancienne energie : " << NRJ_P_coordo << std::endl ;
        std::cout << "Nouvelle energie : " << NRJ_nouveau_vec << std::endl ;
        std::cout.rdbuf(fichier_log_positions.rdbuf());
        affichage_vecteur(P_coordo, n_atomes) ;
        pas/=decremente;
      }

    }
  std::cout.rdbuf(coutbuf);
  fichier_log_NRJ.close();
  fichier_log_positions.close();  
}


void affichage_vecteur(double*& P_coordonnees, const int & n_atomes)
{
  for(int i = 0 ; i < n_atomes ; i++)
    {
      std::cout << std::setw(15) << P_coordonnees[i*2] << " " << P_coordonnees[i*2 + 1] << " ";
    }
  std::cout << std::endl ;
}
