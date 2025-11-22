//////////////////////////////////////
//MINIMISEUR AGREGATS - MINIMISATION//
//////////////////////////////////////

#include "header.h"

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

