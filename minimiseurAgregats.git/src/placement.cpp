/////////////////////////////////////
//FONCTIONS DE PLACEMENT DES ATOMES//
/////////////////////////////////////

#include "header.h"


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
  int n = n_atomes;
  double r = param_supp;
  double a = 0.0;
  double d = 0.0;
  double mT = 0.0;
  double dTheta = 0.0;
  double dPhi = 0.0;

  std::vector <double> vecInter;

  computeValuesFromDatas(n, r, a, d, mT, dTheta, dPhi);

  vecInter = computeResults(r, mT, dPhi);

  for(int i=0; i<vecInter.size()/3; i++)
    {
      P_coordonnees[3*i]=vecInter[3*i];
      P_coordonnees[1+3*i]=vecInter[1+3*i];
      P_coordonnees[2+3*i]=vecInter[2+3*i];
    }
}

void computeValuesFromDatas(const int& n, const double& r, double& a, double& d, double& mT, double& dTheta, double& dPhi)
{
  a = (4.0*M_PI)/(1.0*n) ;
  d = sqrt(a);
  mT = round(M_PI/d) ;
  dTheta = M_PI/mT;
  dPhi = a/dTheta;
}

double placeX(double& theta, double& phi, const double& r)
{
  double results;
  results = sin(theta)*cos(phi);
  return results;
}

double placeY(double& theta, double& phi, const double& r)
{
  double results;
  results = sin(theta)*sin(phi);
  return results;
}

double placeZ(double& theta, const double& r)
{
  double results;
  results = cos(theta);
  return results;
}

std::vector <double> computeResults(double& r, double& mTheta, double& dPhi)
{
  std::vector <double> results;
  double mPhi = 0.0;
  double theta = 0.0;
  double phi = 0.0;
  double r_i = 1.0;

  for(int i=0; i<rint(mTheta); i++)
    {
      theta = (M_PI*(i+0.5))/mTheta;
      mPhi = round((2.0*M_PI*sin(theta))/dPhi);

      for(int j=0; j < rint(mPhi); j++)
        {
          phi = (2.0*M_PI*j)/mPhi;
          results.push_back(placeX(theta, phi, r_i));
          results.push_back(placeY(theta, phi, r_i));
          results.push_back(placeZ(theta, r_i));
        }
      }

  if(r!=1.0)
    {
      for(int i=0; i<results.size()/3; i++)
        {
          results[i*3] *=r;
          results[1+i*3] *=r;
          results[2+i*3] *=r;
        }
    }

  return results;
}


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

  else if(type_forme == "Sphere")
  {
    placer_sphere(n_atomes, param_supp, P_coordonnees);
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