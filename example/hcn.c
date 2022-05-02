#include "lime.h"

double beta= 1.042e-5;
double betahcn= 1.5e-5;
double rnuc = 2.5e2;
double abund = 0.001;

double openAngle = 30. * PI/180.;

double tkin1 = 50.;
double tkin2 = 50.;
double Q1 = 2e28;
double Q2 = 1e28;
double vexp1 = 700.;
double vexp2 = 500.;

/******************************************************************************/

void
input(inputPars *par, image *img){
/*
 * Basic parameters. See cheat sheet for details.
 */
  par->Q1 = Q1;
  par->Q2 = Q2;
  par->openAngle = openAngle;
  par->xne = 1.0;
  par->rHelio       = 1.0;
  par->radius           = 2e8;
  par->minScale         = rnuc;
  par->pIntensity = 1000;
  par->sinkPoints = 500;
  par->moldatfile[0]    = "hcn.dat";
  par->girdatfile[0]    = "g_hcn_1au.dat";
  par->girScale = 1.0;
  par->lte_only         = 0;
  par->useEP       = 0;
  par->nSolveIters  = 1;
  par->outputfile = "hcn.pop";
  par->gridfile         = "grid.vtu";
  par->collPartIds[0]   = 1;
  par->nMolWeights[0]   = 1.0;

/*
 * Definitions for image #0. Add blocks for additional images.
 */
  img[0].velres         = 100.;   // Channel resolution in m/s
  img[0].nchan          = 50;     // Number of channels
  img[0].trans          = 3;    // zero-indexed J quantum number
  img[0].pxls           = 256;    // Pixels per dimension
  img[0].imgres         = 0.5;    // Resolution in arc seconds
  img[0].distance = 1.0*AU; // source distance in m
  img[0].unit           = 0;    // 0:Kelvin 1:Jansky/pixel 2:SI 3:Lsun/pixel 4:tau
  img[0].filename = "hcn.fits"; // Output filename
}

/******************************************************************************/

void
density(double x, double y, double z, double *density){
/*
 * Define variable for radial coordinate
 */
  double r,b,angle;

  b = sqrt(x*x+y*y);
  angle = atan2(b,-z);

  const double rMin = rnuc; /* This cutoff should be chosen smaller than par->minScale but greater than zero (to avoid a singularity at the origin). */

  /*
   * Calculate radial distance from origin
   */
  r=sqrt(x*x+y*y+z*z);
  /*
   * Calculate a Haser density profile
   * (Multiply with 1e6 to go to SI-units)
   */

  if(r<rMin)
      density[0] = 1e-20; /* Just to prevent overflows at r==0! */
  else if(angle<openAngle)
    density[0] = Q1 /(4*PI*pow(r, 2)*vexp1)*exp(-r*beta/vexp1);
  else
    density[0] = Q2 /(4*PI*pow(r, 2)*vexp2)*exp(-r*beta/vexp2);
  
}

/******************************************************************************/

void
temperature(double x, double y, double z, double *temperature){

  double b,angle;

  b = sqrt(x*x+y*y);
  angle = atan2(b,-z);

  if(angle<openAngle){
    temperature[0] = tkin1;
  }
  else
    temperature[0] = tkin2;
}

/******************************************************************************/

void
molNumDensity(double x, double y, double z, double *nmol){
 /*
 * Define variable for radial coordinate
 */
  double r,b,angle;

  b = sqrt(x*x+y*y);
  angle = atan2(b,-z);

  const double rMin = rnuc; /* This cutoff should be chosen smaller than par->minScale but greater than zero (to avoid a singularity at the origin). */

  /*
   * Calculate radial distance from origin
   */
  r=sqrt(x*x+y*y+z*z);
  /*
   * Calculate a Haser density profile
   * (Multiply with 1e6 to go to SI-units)
   */

  if(r<rMin)
    nmol[0] = 0.;
  else if(angle<openAngle)
    nmol[0] =abund*Q1/(4*PI*pow(r, 2)*vexp1)*exp(-r*betahcn/vexp1);
  else
    nmol[0] =abund*Q2/(4*PI*pow(r, 2)*vexp2)*exp(-r*betahcn/vexp2);
}

/******************************************************************************/

void
doppler(double x, double y, double z, double *doppler){
  *doppler = 100.;
}

/******************************************************************************/

void
velocity(double x, double y, double z, double *vel){
/*
 * Variables for spherical coordinates
 */
  double phi, theta,b,angle;

  b = sqrt(x*x+y*y);
  angle = atan2(b,-z);
/*
 * Transform Cartesian coordinates into spherical coordinates
 */
  theta=atan2(sqrt(x*x+y*y),z);
  phi=atan2(y,x);
/*
 * Vector transformation back into Cartesian basis
 */
  if(angle<openAngle){
    vel[0]=vexp1*sin(theta)*cos(phi);
    vel[1]=vexp1*sin(theta)*sin(phi);
    vel[2]=vexp1*cos(theta);

  }
  else{
    vel[0]=vexp2*sin(theta)*cos(phi);
    vel[1]=vexp2*sin(theta)*sin(phi);
    vel[2]=vexp2*cos(theta);
  }
}

