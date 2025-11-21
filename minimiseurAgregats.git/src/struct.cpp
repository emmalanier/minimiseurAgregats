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