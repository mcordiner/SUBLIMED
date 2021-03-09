#include "lime.h"

double beta= 1.042e-5;
double betahcn= 1.5e-5;

double tkin = 50.;
double rnuc = 2.5e2;
double Q1 = 2e28;
double Q2 = 1e28;
double vexp1 = 700.;
double vexp2 = 500.;

double abund = 0.001;

/******************************************************************************/

void
input(inputPars *par, image *img){
/*
 * Basic parameters. See cheat sheet for details.
 */
  par->Q1 = Q1;
  par->Q2 = Q2;
  par->openAng = 30. * PI/180.
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
  par->gridfile         = "grid.vtk";
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
density1(double x, double y, double z, double *density){

  double r;

  r=sqrt(x*x+y*y+z*z);

  if(r<rnuc)
    density[0] = 1e-20; /* Just to prevent overflows at r==0! */
  else
    density[0] =  /(4.*PI*pow(r, 2)*vexp)*exp(-r*beta/vexp);
}

void
density2(double x, double y, double z, double *density){

  double r;

  r=sqrt(x*x+y*y+z*z);

  if(r<rnuc)
    density[0] = 1e-20; /* Just to prevent overflows at r==0! */
  else
    density[0] =  /(4.*PI*pow(r, 2)*vexp)*exp(-r*beta/vexp);
}

/******************************************************************************/

void
temperature1(double x, double y, double z, double *temperature){
  temperature[0] = tkin;
}

void
temperature2(double x, double y, double z, double *temperature){
  temperature[0] = tkin;
}

/******************************************************************************/

void
molNumDensity1(double x, double y, double z, double *nmol){
  
  double r;
  
  r=sqrt(x*x+y*y+z*z);

  if(r<rnuc)
    nmol[0] = 0.;
  else
    nmol[0] = abund*Q1/(4*PI*pow(r, 2)*vexp1)*exp(-r*betahcn/vexp1);
}

void
molNumDensity2(double x, double y, double z, double *nmol){
  
  double r;
  
  r=sqrt(x*x+y*y+z*z);

  if(r<rnuc)
    nmol[0] = 0.;
  else
    nmol[0] = abund*Q2/(4*PI*pow(r, 2)*vexp2)*exp(-r*betahcn/vexp2);
}


/******************************************************************************/

void
doppler(double x, double y, double z, double *doppler){
  *doppler = 100.;
}

/******************************************************************************/

void
velocity1(double x, double y, double z, double *vel){
/*
 * Variables for spherical coordinates
 */
  double phi, theta;
/*
 * Transform Cartesian coordinates into spherical coordinates
 */
  theta=atan2(sqrt(x*x+y*y),z);
  phi=atan2(y,x);
/*
 * Vector transformation back into Cartesian basis
 */
  vel[0]=vexp1*sin(theta)*cos(phi);
  vel[1]=vexp1*sin(theta)*sin(phi);
  vel[2]=vexp1*cos(theta);
}

void
velocity2(double x, double y, double z, double *vel){
/*
 * Variables for spherical coordinates
 */
  double phi, theta;
/*
 * Transform Cartesian coordinates into spherical coordinates
 */
  theta=atan2(sqrt(x*x+y*y),z);
  phi=atan2(y,x);
/*
 * Vector transformation back into Cartesian basis
 */
  vel[0]=vexp2*sin(theta)*cos(phi);
  vel[1]=vexp2*sin(theta)*sin(phi);
  vel[2]=vexp2*cos(theta);
}

