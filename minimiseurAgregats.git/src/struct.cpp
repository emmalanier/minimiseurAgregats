////////////////////////////////////////////////////////////////
//MINIMISEUR AGREGATS - SURCHARGES ET METHODES POUR STRUCTURES//
////////////////////////////////////////////////////////////////

#include "header.h"

vecteur operator-(const vecteur& va, const vecteur& vb)
{
  vecteur results;

  results.compoX = va.compoX - vb.compoX;
  results.compoY = va.compoY - vb.compoY;
  results.compoZ = va.compoZ - vb.compoZ;

  return results;
}

vecteur operator-(const position& pa, const position& pb)
{
  vecteur results;

  results.compoX = pa.x - pb.x;
  results.compoY = pa.y - pb.y;
  results.compoZ = pa.z - pb.z;

  return results;
}

vecteur operator+(const vecteur& va, const vecteur& vb)
{
  vecteur results;
  
  results.compoX = va.compoX + vb.compoX;
  results.compoY = va.compoY + vb.compoY;
  results.compoZ = va.compoZ + vb.compoZ;

  return results;
}

double operator*(const vecteur& va, const vecteur& vb)
{
  double results;

  results = va.compoX*vb.compoX + va.compoY*vb.compoY + va.compoZ*vb.compoZ;

  return results;
}

vecteur operator*(const double& scalaire, const vecteur& va)
{
  vecteur results;

  results.compoX = va.compoX * scalaire;
  results.compoY = va.compoY * scalaire;
  results.compoZ = va.compoZ * scalaire;

  return results;
}

position operator+(const position& p, const vecteur& va)
{
  position results ;

  results.x = p.x + va.compoX;
  results.y = p.y + va.compoY;
  results.z = p.z + va.compoZ;

  return results ;
}


vecteur vecteur::normalize()
{
  vecteur results;

  double norme = sqrt(compoX*compoX + compoY*compoY + compoZ*compoZ);

  results.compoX = compoX/norme;
  results.compoY = compoY/norme;
  results.compoZ = compoZ/norme;

  return results;
}

void force::setToZero()
{
  valeur = 0.0;
  
  vecForce.compoX = 0.0;
  vecForce.compoY = 0.0;
  vecForce.compoZ = 0.0;
  
  vecForce.ptApp.x = 0.0;
  vecForce.ptApp.y = 0.0;
  vecForce.ptApp.z = 0.0;
}

void force::operator+=(const force& f)
{
  valeur += f.valeur;
  vecForce = vecForce + f.vecForce;
  
  vecForce.ptApp.x = vecForce.ptApp.x + f.vecForce.ptApp.x ;
  vecForce.ptApp.y = vecForce.ptApp.y + f.vecForce.ptApp.y ;
  vecForce.ptApp.z = vecForce.ptApp.z + f.vecForce.ptApp.z ;
}

force operator*(const double& scalaire, const force& f)
{
  force results;

  results.vecForce.compoX = results.vecForce.compoX * scalaire;
  results.vecForce.compoY = results.vecForce.compoY * scalaire;
  results.vecForce.compoZ = results.vecForce.compoZ * scalaire;

  return results;
}