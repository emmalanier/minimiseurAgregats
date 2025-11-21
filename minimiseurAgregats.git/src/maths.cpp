///////////////////////////////
//MINIMISEUR AGREGATS - MATHS//
///////////////////////////////

#include "fonctions.h"

////////////////////
//CALCUL VECTORIEL//
////////////////////

double calcul_distance_3D(double x, double y, double z)
{
  double results = 0.0;
  double results_inter = 0.0;

  results_inter = (x*x)+(y*y)+(z*z);

  results = sqrt(results_inter);

  return results;
}

//////////////
//POTENTIELS//
//////////////

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

///////////
//ENERGIE//
///////////

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

//////////
//FORCES//
//////////

force calcForceElec(partChargee part1, partChargee part2)
{
  force results;

  vecteur r = part2.lieu - part1.lieu;
  r.norme = sqrt(r.compoX*r.compoX + r.compoY*r.compoY + r.compoZ*r.compoZ;
  double r2 = r.norme*r.norme;
  vecteur direction = distance1vers2.normalize();

  results.valeur = (part1.charge*part2.charge)/(4*M_PI*epsilon_0*distance1Vers2Carre);

  results.compoX = results.valeur * direction.compoX;
  results.compoY = results.valeur * direction.compoY;
  results.compoZ = results.valeur * direction.compoZ;

  return results;

}

vecteur calculAcceleration(std::vector <force> & vecForces, partChargee & particule)
{
  //On part du principe que nous sommes dans une cadre non-relativiste, donc F = ma

  vecteur results;

  vecteur sommeForces;

  sommeForces.compoX = 0.0;
  sommeForces.compoY = 0.0;
  sommeForces.compoZ = 0.0;

  for(int i=0; i<vecForces.size(); i++)
    sommeForces += vecForces[i];

  results = sommeForces*(1.0/particule.masse);

  return results;
}
  
void partChargee::update(double dt);
{
  vitesse = vitesse + acceleration*dt;
  lieu = lieu + vitesse*dt;
}
