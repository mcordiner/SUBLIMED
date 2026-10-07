/*
 *  solver.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015-2017 The LIME development team
 *
TODO:
  - The test to run calculateJBar() etc in levelPops just tests dens[0]. This is a bit sloppy.
 */

#include <stdio.h>
#include <stdlib.h>

#include <cvode/cvode.h>               /* prototypes for CVODE fcts., consts.  */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sunlinsol/sunlinsol_spgmr.h>   /* access to SPGMR SUNLinearSolver       */
#include <sunmatrix/sunmatrix_dense.h>   /* access to dense SUNMatrix             */
#include <sunlinsol/sunlinsol_dense.h>   /* access to dense SUNLinearSolver       */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */

#include "lime.h"
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_bessel.h>

//###
#define Ith(v,i)    NV_Ith_S(v,i)         /* Ith numbers components 0..NEQ-1 */
#define IJth(sunMatrix,i,j) SM_ELEMENT_D(sunMatrix,i,j) /* IJth numbers rows,cols 0..NEQ-1 */

/* Data concerning a single grid vertex which is passed from calculateJBar() to solveStatEq(). This data needs to be thread-safe. */
typedef struct {
  double *jbar,*phot,*vfac,*vfac_loc;
} gridPointData;

struct blend{
  int molJ, lineJ;
  double deltaV;
};

struct lineWithBlends{
  int lineI, numBlends;
  struct blend *blends;
};

struct molWithBlends{
  int molI, numLinesWithBlends;
  struct lineWithBlends *lines;
};

struct blendInfo{
  int numMolsWithBlends;
  struct molWithBlends *mols;
};
struct time_struct{
  int *id;
  double *time;
};


/* Parameters used to determine transition rates */
struct transitionParams{
  int array_size; //Equal to the number of levels (i.e NEQ)
  double *A_array; //Holds Einstein's As
  molData *md;
  int ispec;
  struct grid *gp;
  configInfo *par;
  double *jbar_grid;
  int *nMaserWarnings;
  int *gp_sorter;
  int subGrid;
  int gp_pIntensity;
  
  /* --- NEW CACHING & MEMORY VARIABLES --- */
  double *p_rates;      // Working matrix allocated once
  double *coll_rates;   // Radius-dependent rates (cached)
  double last_radius;   // Track when radius changes
  double vexp;          // Cached expansion velocity
};

/* Checks for errors when calling any CVode functions */
int check_retval(void *returnvalue, const char *funcname, int opt)
{
  int *retval;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && returnvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
      funcname);
    return(1); }

  /* Check if retval < 0 */
  else if (opt == 1) {
    retval = (int *) returnvalue;
    if (*retval < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
        funcname, *retval);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && returnvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
      funcname);
    return(1); }

  return(0);
}


/*ELECTRON TEMPERATURE FUNCTION*/
double Telec(double r, double Q, double Tkin){
  double Te,rcs;
  double Tmax = 1e4;
  rcs = 1.125e6 * pow(Q/1e29,0.75);
  if (r < rcs) {
    Te = Tkin;
  }
  else if (r > 2.*rcs) {
    Te = Tmax;
  }
  else {
    Te = Tkin + (Tmax - Tkin)*((r/rcs)-1.);
  }
  return Te;
}

/*ELECTRON DENSITY FUNCTION*/
double nelec(double r, double Q, double vexp, double Te, double rH, double xne){
  double Rrec, ne, krec, kion;
  kion = 4.1e-7;
  krec = 3e-13 * sqrt(300./Te); /* Recombination rate From RATE12, accountiong for cm3 to m3 conversion */
  Rrec = 3.2e6 * sqrt(Q/1e29);
  /* Equation 5 of Zakharov 2007 */
  ne = xne * sqrt(Q*kion/vexp/krec/(rH*rH)) * pow((Te/300.),0.15) * (Rrec/(r*r)) * (1-exp(-r/Rrec)) + (5e6/(rH*rH));
  
  return ne;
}



/*....................................................................*/
void
freeMolsWithBlends(struct molWithBlends *mols, const int numMolsWithBlends){
  int mi, li;

  if(mols != NULL){
    for(mi=0;mi<numMolsWithBlends;mi++){
      if(mols[mi].lines != NULL){
        for(li=0;li<mols[mi].numLinesWithBlends;li++)
          free(mols[mi].lines[li].blends);
        free(mols[mi].lines);
      }
    }
    free(mols);
  }
}

/*....................................................................*/
void
freeGridPointData(const int nSpecies, gridPointData *mol){
  /*
Note that this is called from within the multi-threaded block.
  */
  int i;
  if(mol!= NULL){
    for(i=0;i<nSpecies;i++){
      free(mol[i].jbar);
      free(mol[i].phot);
      free(mol[i].vfac);
      free(mol[i].vfac_loc);
    }
    free(mol);
  }
}

/*....................................................................*/
void lineBlend(molData *m, configInfo *par, struct blendInfo *blends){
  /*
This obtains information on all the lines of all the radiating species which have other lines within some cutoff velocity separation.
  */
  int molI, lineI, molJ, lineJ;
  int nmwb, nlwb, numBlendsFound, li, bi;
  double deltaV;
  struct blend *tempBlends=NULL;
  struct lineWithBlends *tempLines=NULL;

  (*blends).mols = malloc(sizeof(struct molWithBlends)*par->nSpecies);
  (*blends).numMolsWithBlends = 0;

  nmwb = 0;
  for(molI=0;molI<par->nSpecies;molI++){
    tempBlends = malloc(sizeof(struct blend)*m[molI].nline);
    tempLines  = malloc(sizeof(struct lineWithBlends)*m[molI].nline);

    nlwb = 0;
    for(lineI=0;lineI<m[molI].nline;lineI++){
      numBlendsFound = 0;
      for(molJ=0;molJ<par->nSpecies;molJ++){
        for(lineJ=0;lineJ<m[molJ].nline;lineJ++){
          if(!(molI==molJ && lineI==lineJ)){
            deltaV = (m[molJ].freq[lineJ] - m[molI].freq[lineI])*CLIGHT/m[molI].freq[lineI];
            if(fabs(deltaV)<maxBlendDeltaV){
              tempBlends[numBlendsFound].molJ   = molJ;
              tempBlends[numBlendsFound].lineJ  = lineJ;
              tempBlends[numBlendsFound].deltaV = deltaV;
              numBlendsFound++;
            }
          }
        }
      }

      if(numBlendsFound>0){
        tempLines[nlwb].lineI = lineI;
        tempLines[nlwb].numBlends = numBlendsFound;
        tempLines[nlwb].blends = malloc(sizeof(struct blend)*numBlendsFound);
        for(bi=0;bi<numBlendsFound;bi++)
          tempLines[nlwb].blends[bi] = tempBlends[bi];

        nlwb++;
      }
    }

    if(nlwb>0){
      (*blends).mols[nmwb].molI = molI;
      (*blends).mols[nmwb].numLinesWithBlends = nlwb;
      (*blends).mols[nmwb].lines = malloc(sizeof(struct lineWithBlends)*nlwb);
      for(li=0;li<nlwb;li++){
        (*blends).mols[nmwb].lines[li].lineI     = tempLines[li].lineI;
        (*blends).mols[nmwb].lines[li].numBlends = tempLines[li].numBlends;
        (*blends).mols[nmwb].lines[li].blends = malloc(sizeof(struct blend)*tempLines[li].numBlends);
        for(bi=0;bi<tempLines[li].numBlends;bi++)
          (*blends).mols[nmwb].lines[li].blends[bi] = tempLines[li].blends[bi];
      }

      nmwb++;
    }

    free(tempLines);
    free(tempBlends);
  }

  (*blends).numMolsWithBlends = nmwb;
  if(nmwb>0){
    if(!par->blend)
      if(!silent) warning("Blended lines are present in the model.");

    (*blends).mols = realloc((*blends).mols, sizeof(struct molWithBlends)*nmwb);
  }else{
    if(par->blend)
      if(!silent) warning("Line blending is switched on, but no blended lines were found.");

    free((*blends).mols);
    (*blends).mols = NULL;
  }
}

/*....................................................................*/
void calcGridCollRates(configInfo *par, molData *md, struct grid *gp){
  int i,id,ipart,itrans,itemp,tnint=-1;
  struct cpData part;
  double fac;

  for(i=0;i<par->nSpecies;i++){
    for(id=0;id<par->ncell;id++){
      gp[id].mol[i].partner = malloc(sizeof(struct rates)*md[i].npart);
    }

    for(ipart=0;ipart<md[i].npart;ipart++){
      part = md[i].part[ipart];
      for(id=0;id<par->ncell;id++){
        for(itrans=0;itrans<part.ntrans;itrans++){
          if((gp[id].t[0]>part.temp[0])&&(gp[id].t[0]<part.temp[part.ntemp-1])){
            for(itemp=0;itemp<part.ntemp-1;itemp++){
              if((gp[id].t[0]>part.temp[itemp])&&(gp[id].t[0]<=part.temp[itemp+1])){
                tnint=itemp;
              }
            }
            fac=(gp[id].t[0]-part.temp[tnint])/(part.temp[tnint+1]-part.temp[tnint]);
            gp[id].mol[i].partner[ipart].t_binlow = tnint;
            gp[id].mol[i].partner[ipart].interp_coeff = fac;

    } else if(gp[id].t[0]<=part.temp[0]) {
      gp[id].mol[i].partner[ipart].t_binlow = 0;
      gp[id].mol[i].partner[ipart].interp_coeff = 0.0;
    } else {
      gp[id].mol[i].partner[ipart].t_binlow = part.ntemp-2;
      gp[id].mol[i].partner[ipart].interp_coeff = 1.0;
    }
        } /* End loop over transitions. */
      } /* End loop over grid points. */
    } /* End loop over collision partners. */
  } /* End loop over radiating molecules. */
}

/*....................................................................*/
void mallocGridCont(configInfo *par, molData *md, struct grid *gp){
  int id,si,li;

  for(id=0;id<par->ncell;id++){
    for(si=0;si<par->nSpecies;si++){
      gp[id].mol[si].cont = malloc(sizeof(*(gp[id].mol[si].cont))*md[si].nline);
      for(li=0;li<md[si].nline;li++){
        gp[id].mol[si].cont[li].dust = 0.0;
        gp[id].mol[si].cont[li].knu  = 0.0;
      }
    }
  }
}

/*....................................................................*/
void freeGridCont(configInfo *par, struct grid *gp){
  int id,si;

  for(id=0;id<par->ncell;id++){
    if(gp[id].mol==NULL)
      continue;

    for(si=0;si<par->nSpecies;si++){
      free(gp[id].mol[si].cont);
      gp[id].mol[si].cont = NULL;
    }
  }
}

/*....................................................................*/
void calcGridLinesDustOpacity(configInfo *par, molData *md, double *lamtab\
  , double *kaptab, const int nEntries, struct grid *gp){

  int iline,id,si;
  double *kappatab,gtd;
  gsl_spline *spline = NULL;
  gsl_interp_accel *acc = NULL;
  double *knus=NULL, *dusts=NULL;

  if(par->dust != NULL){
    acc = gsl_interp_accel_alloc();
    spline = gsl_spline_alloc(gsl_interp_cspline,nEntries);
    gsl_spline_init(spline,lamtab,kaptab,nEntries);
  }

  for(si=0;si<par->nSpecies;si++){
    kappatab = malloc(sizeof(*kappatab)*md[si].nline);
    knus     = malloc(sizeof(*knus)    *md[si].nline);
    dusts    = malloc(sizeof(*dusts)   *md[si].nline);

    if(par->dust == NULL){
      for(iline=0;iline<md[si].nline;iline++)
        kappatab[iline] = 0.;
    }else{
      for(iline=0;iline<md[si].nline;iline++)
        kappatab[iline] = interpolateKappa(md[si].freq[iline]\
                        , lamtab, kaptab, nEntries, spline, acc);
    }

    for(id=0;id<par->ncell;id++){
      gasIIdust(gp[id].x[0],gp[id].x[1],gp[id].x[2],&gtd);
      calcDustData(par, gp[id].dens, md[si].freq, gtd, kappatab, md[si].nline, gp[id].t, knus, dusts);
      for(iline=0;iline<md[si].nline;iline++){
        gp[id].mol[si].cont[iline].knu  = knus[iline];
        gp[id].mol[si].cont[iline].dust = dusts[iline];
      }
    }

    free(kappatab);
    free(knus);
    free(dusts);
  }

  if(par->dust != NULL){
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
  }
}

/*....................................................................*/
int
getNextEdge(double *inidir, const int startGi, const int presentGi\
  , struct grid *gp, const gsl_rng *ran){
  /*
Note that this is called from within the multi-threaded block.
  */
  int i,ni,niOfSmallest=-1,niOfNextSmallest=-1;
  double dirCos,distAlongTrack,dirFromStart[3],coord,distToTrackSquared,smallest=0.0,nextSmallest=0.0;
  const static double scatterReduction = 0.4;

  i = 0;
  for(ni=0;ni<gp[presentGi].numNeigh;ni++){
    dirCos = dotProduct3D(inidir, gp[presentGi].dir[ni].xn);

    if(dirCos<=0.0)
  continue; /* because the edge points in the backward direction. */

    dirFromStart[0] = gp[presentGi].neigh[ni]->x[0] - gp[startGi].x[0];
    dirFromStart[1] = gp[presentGi].neigh[ni]->x[1] - gp[startGi].x[1];
    dirFromStart[2] = gp[presentGi].neigh[ni]->x[2] - gp[startGi].x[2];
    distAlongTrack = dotProduct3D(inidir, dirFromStart);

    coord = dirFromStart[0] - distAlongTrack*inidir[0];
    distToTrackSquared  = coord*coord;
    coord = dirFromStart[1] - distAlongTrack*inidir[1];
    distToTrackSquared += coord*coord;
    coord = dirFromStart[2] - distAlongTrack*inidir[2];
    distToTrackSquared += coord*coord;

    if(i==0){
      smallest = distToTrackSquared;
      niOfSmallest = ni;
    }else{
      if(distToTrackSquared<smallest){
        nextSmallest = smallest;
        niOfNextSmallest = niOfSmallest;
        smallest = distToTrackSquared;
        niOfSmallest = ni;
      }else if(i==1 || distToTrackSquared<nextSmallest){
        nextSmallest = distToTrackSquared;
        niOfNextSmallest = ni;
      }
    }

    i++;
  }

  if(i>1){
    if((smallest + scatterReduction*nextSmallest)*gsl_rng_uniform(ran)<smallest){
      return niOfNextSmallest;
    }else{
      return niOfSmallest;
    }
  }else if(i>0){
    return niOfSmallest;
  }else{
    if(!silent)
      bail_out("Photon propagation error - no valid edges.");
    exit(1);
  }
}

/*....................................................................*/
void calcLineAmpPWLin(struct grid *g, const int id, const int k\
  , const int molI, const double deltav, double *inidir, double *vfac_in, double *vfac_out){
  double binv_this, binv_next, v[5];

  binv_this=g[id].mol[molI].binv;
  binv_next=(g[id].neigh[k])->mol[molI].binv;
  v[0]=deltav-dotProduct3D(inidir,g[id].vel);
  v[1]=deltav-dotProduct3D(inidir,&(g[id].v1[3*k]));
  v[2]=deltav-dotProduct3D(inidir,&(g[id].v2[3*k]));
  v[3]=deltav-dotProduct3D(inidir,&(g[id].v3[3*k]));
  v[4]=deltav-dotProduct3D(inidir,g[id].neigh[k]->vel);

  if (fabs(v[1]-v[0])*binv_this>(2.0*BIN_WIDTH)) {
     *vfac_out=0.5*geterf(v[0]*binv_this,v[1]*binv_this);
  } else *vfac_out=0.5*gaussline(0.5*(v[0]+v[1]),binv_this);
  if (fabs(v[2]-v[1])*binv_this>(2.0*BIN_WIDTH)) {
    *vfac_out+=0.5*geterf(v[1]*binv_this,v[2]*binv_this);
  } else *vfac_out+=0.5*gaussline(0.5*(v[1]+v[2]),binv_this);

  if (fabs(v[3]-v[2])*binv_next>(2.0*BIN_WIDTH)) {
     *vfac_in=0.5*geterf(v[2]*binv_next,v[3]*binv_next);
  } else *vfac_in=0.5*gaussline(0.5*(v[2]+v[3]),binv_next);
  if (fabs(v[4]-v[3])*binv_next>(2.0*BIN_WIDTH)) {
    *vfac_in+=0.5*geterf(v[3]*binv_next,v[4]*binv_next);
  } else *vfac_in+=0.5*gaussline(0.5*(v[3]+v[4]),binv_next);
}

/*....................................................................*/
void calcLineAmpLin(struct grid *g, const int id, const int k\
  , const int molI, const double deltav, double *inidir, double *vfac_in, double *vfac_out){
  double binv_this, binv_next, v[3];

  binv_this=g[id].mol[molI].binv;
  binv_next=(g[id].neigh[k])->mol[molI].binv;
  v[0]=deltav-dotProduct3D(inidir,g[id].vel);
  v[2]=deltav-dotProduct3D(inidir,g[id].neigh[k]->vel);
  v[1]=0.5*(v[0]+v[2]);

  if (fabs(v[1]-v[0])*binv_this>(2.0*BIN_WIDTH)) {
     *vfac_out=geterf(v[0]*binv_this,v[1]*binv_this);
  } else *vfac_out+=gaussline(0.5*(v[0]+v[1]),binv_this);

  if (fabs(v[2]-v[1])*binv_next>(2.0*BIN_WIDTH)) {
     *vfac_in=geterf(v[1]*binv_next,v[2]*binv_next);
  } else *vfac_in+=gaussline(0.5*(v[1]+v[2]),binv_next);
}

/*....................................................................*/
void
calculateJBar(int id, struct grid *gp, molData *md, const gsl_rng *ran\
  , configInfo *par, const int nlinetot, struct blendInfo blends\
  , gridPointData *mp, double *halfFirstDs, int *nMaserWarnings){
  /*
Note that this is called from within the multi-threaded block.
  */

  int iphot,iline,here,there,firststep,neighI,numLinks=0;
  int nextMolWithBlend, nextLineWithBlend, molI, lineI, molJ, lineJ, bi;
  double segment,vblend_in,vblend_out,dtau,expDTau,ds_in=0.0,ds_out=0.0,pt_theta,pt_z,semiradius;
  double deltav[par->nSpecies],vfac_in[par->nSpecies],vfac_out[par->nSpecies],vfac_inprev[par->nSpecies];
  double expTau[nlinetot],inidir[3];
  double remnantSnu,velProj;
  char message[STR_LEN_0];

  for(iphot=0;iphot<gp[id].nphot;iphot++){
    firststep=1;
    iline = 0;
    for(molI=0;molI<par->nSpecies;molI++){
      for(lineI=0;lineI<md[molI].nline;lineI++){
        mp[molI].phot[lineI+iphot*md[molI].nline]=0.;
        expTau[iline]=1.;
        iline++;
      }
    }

    pt_theta=gsl_rng_uniform(ran)*2*M_PI;
    pt_z=2*gsl_rng_uniform(ran)-1;
    semiradius = sqrt(1.-pt_z*pt_z);
    inidir[0]=semiradius*cos(pt_theta);
    inidir[1]=semiradius*sin(pt_theta);
    inidir[2]=pt_z;

    segment=gsl_rng_uniform(ran)-0.5;

    for (molI=0;molI<par->nSpecies;molI++){
      deltav[molI]=4.3*segment*gp[id].mol[molI].dopb+dotProduct3D(inidir,gp[id].vel);
      mp[molI].vfac_loc[iphot]=gaussline(deltav[molI]-dotProduct3D(inidir,gp[id].vel),gp[id].mol[molI].binv);
    }

    here = gp[id].id;

    numLinks=0;
    while(!gp[here].sink){
      numLinks++;
      if(numLinks>par->ncell){
        if(!silent){
          snprintf(message, STR_LEN_0, "Bad grid? Too many links in photon path, point %d photon %d", id, iphot);
          bail_out(message);
        }
exit(1);
      }

      neighI = getNextEdge(inidir,id,here,gp,ran);

      there=gp[here].neigh[neighI]->id;

      if(firststep){
        firststep=0;
        ds_out=0.5*gp[here].ds[neighI]*dotProduct3D(inidir,gp[here].dir[neighI].xn);
        halfFirstDs[iphot]=ds_out;

        for(molI=0;molI<par->nSpecies;molI++){
          if(par->edgeVelsAvailable) {
            calcLineAmpPWLin(gp,here,neighI,molI,deltav[molI],inidir,&vfac_in[molI],&vfac_out[molI]);
         } else
            calcLineAmpLin(gp,here,neighI,molI,deltav[molI],inidir,&vfac_in[molI],&vfac_out[molI]);

          mp[molI].vfac[iphot]=vfac_out[molI];
        }
        here=there;
    continue;
      }

      ds_in=ds_out;
      ds_out=0.5*gp[here].ds[neighI]*dotProduct3D(inidir,gp[here].dir[neighI].xn);

      for(molI=0;molI<par->nSpecies;molI++){
        vfac_inprev[molI]=vfac_in[molI];
        if(par->edgeVelsAvailable)
          calcLineAmpPWLin(gp,here,neighI,molI,deltav[molI],inidir,&vfac_in[molI],&vfac_out[molI]);
        else
          calcLineAmpLin(gp,here,neighI,molI,deltav[molI],inidir,&vfac_in[molI],&vfac_out[molI]);
      }

      nextMolWithBlend = 0;
      iline = 0;
      for(molI=0;molI<par->nSpecies;molI++){
        nextLineWithBlend = 0;
        for(lineI=0;lineI<md[molI].nline;lineI++){
          double jnu_line_in=0., jnu_line_out=0., jnu_cont=0., jnu_blend=0.;
          double alpha_line_in=0., alpha_line_out=0., alpha_cont=0., alpha_blend=0.;

          sourceFunc_line(&md[molI],vfac_inprev[molI],&(gp[here].mol[molI]),lineI,&jnu_line_in,&alpha_line_in);
          sourceFunc_line(&md[molI],vfac_out[molI],&(gp[here].mol[molI]),lineI,&jnu_line_out,&alpha_line_out);
          sourceFunc_cont(gp[here].mol[molI].cont[lineI],&jnu_cont,&alpha_cont);

          if(par->blend && blends.mols!=NULL && molI==blends.mols[nextMolWithBlend].molI\
          && lineI==blends.mols[nextMolWithBlend].lines[nextLineWithBlend].lineI){

            for(bi=0;bi<blends.mols[nextMolWithBlend].lines[nextLineWithBlend].numBlends;bi++){
              molJ  = blends.mols[nextMolWithBlend].lines[nextLineWithBlend].blends[bi].molJ;
              lineJ = blends.mols[nextMolWithBlend].lines[nextLineWithBlend].blends[bi].lineJ;
              velProj = deltav[molI] - blends.mols[nextMolWithBlend].lines[nextLineWithBlend].blends[bi].deltaV;

              if(par->edgeVelsAvailable)
                calcLineAmpPWLin(gp,here,neighI,molJ,velProj,inidir,&vblend_in,&vblend_out);
              else
                calcLineAmpLin(gp,here,neighI,molJ,velProj,inidir,&vblend_in,&vblend_out);

              sourceFunc_line(&md[molJ],vblend_out,&(gp[here].mol[molJ]),lineJ,&jnu_blend,&alpha_blend);
            }

            nextLineWithBlend++;
            if(nextLineWithBlend>=blends.mols[nextMolWithBlend].numLinesWithBlends){
              nextLineWithBlend = 0;
            }
          }

    dtau=(alpha_line_out+alpha_cont+alpha_blend)*ds_out;
          if(dtau < -MAX_NEG_OPT_DEPTH) dtau = -MAX_NEG_OPT_DEPTH;
          calcSourceFn(dtau, par, &remnantSnu, &expDTau);
          remnantSnu *= (jnu_line_out+jnu_cont+jnu_blend)*ds_out;
          mp[molI].phot[lineI+iphot*md[molI].nline]+=expTau[iline]*remnantSnu;
    expTau[iline]*=expDTau;

    dtau=(alpha_line_in+alpha_cont+alpha_blend)*ds_in;
          if(dtau < -MAX_NEG_OPT_DEPTH) dtau = -MAX_NEG_OPT_DEPTH;
          calcSourceFn(dtau, par, &remnantSnu, &expDTau);
          remnantSnu *= (jnu_line_in+jnu_cont+jnu_blend)*ds_in;
          mp[molI].phot[lineI+iphot*md[molI].nline]+=expTau[iline]*remnantSnu;
    expTau[iline]*=expDTau;

          if(expTau[iline] > exp(MAX_NEG_OPT_DEPTH)){
            (*nMaserWarnings)++;
            expTau[iline]=exp(MAX_NEG_OPT_DEPTH);
          }

          iline++;
        }

        if(par->blend && blends.mols!=NULL && molI==blends.mols[nextMolWithBlend].molI)
          nextMolWithBlend++;
      }

      here=there;
    };

    iline = 0;
    for(molI=0;molI<par->nSpecies;molI++){
      for(lineI=0;lineI<md[molI].nline;lineI++){
        mp[molI].phot[lineI+iphot*md[molI].nline]+=expTau[iline]*md[molI].cmb[lineI];
        iline++;
      }
    }
  }
}

/*....................................................................*/
void
updateJBar(int posn, molData *md, struct grid *gp, const int molI\
  , configInfo *par, struct blendInfo blends, int nextMolWithBlend\
  , gridPointData *mp, double *halfFirstDs){
  /*
Note that this is called from within the multi-threaded block.
  */
  int lineI,iphot,bi,molJ,lineJ,nextLineWithBlend;
  double dtau,expDTau,remnantSnu,vsum=0.;
  
  for(lineI=0;lineI<md[molI].nline;lineI++) mp[molI].jbar[lineI]=0.;

  for(iphot=0;iphot<gp[posn].nphot;iphot++){
    if(mp[molI].vfac_loc[iphot]>0){
      nextLineWithBlend = 0;
      for(lineI=0;lineI<md[molI].nline;lineI++){
        double jnu=0.0;
        double alpha=0;

        sourceFunc_line(&md[molI],mp[molI].vfac[iphot],&(gp[posn].mol[molI]),lineI,&jnu,&alpha);
        sourceFunc_cont(gp[posn].mol[molI].cont[lineI],&jnu,&alpha);

        if(par->blend && blends.mols!=NULL && molI==blends.mols[nextMolWithBlend].molI\
        && lineI==blends.mols[nextMolWithBlend].lines[nextLineWithBlend].lineI){
          for(bi=0;bi<blends.mols[nextMolWithBlend].lines[nextLineWithBlend].numBlends;bi++){
            molJ  = blends.mols[nextMolWithBlend].lines[nextLineWithBlend].blends[bi].molJ;
            lineJ = blends.mols[nextMolWithBlend].lines[nextLineWithBlend].blends[bi].lineJ;
            sourceFunc_line(&md[molJ],mp[molI].vfac[iphot],&(gp[posn].mol[molJ]),lineJ,&jnu,&alpha);
          }

          nextLineWithBlend++;
          if(nextLineWithBlend>=blends.mols[nextMolWithBlend].numLinesWithBlends){
            nextLineWithBlend = 0;
          }
        }

        dtau=alpha*halfFirstDs[iphot];
        calcSourceFn(dtau, par, &remnantSnu, &expDTau);
        remnantSnu *= jnu*halfFirstDs[iphot];
        mp[molI].jbar[lineI]+=mp[molI].vfac_loc[iphot]*(expDTau*mp[molI].phot[lineI+iphot*md[molI].nline]+remnantSnu);

      }
      vsum+=mp[molI].vfac_loc[iphot];
    }
  }
  for(lineI=0;lineI<md[molI].nline;lineI++) mp[molI].jbar[lineI] /= vsum;
}

/*....................................................................*/
void lteOnePoint(molData *md, const int ispec, const double temp, double *pops){
  int ilev;
  double sum;

  sum = 0.0;
  for(ilev=0;ilev<md[ispec].nlev;ilev++){
    pops[ilev] = md[ispec].gstat[ilev]*exp(-HCKB*md[ispec].eterm[ilev]/temp);
    sum += pops[ilev];
  }
  for(ilev=0;ilev<md[ispec].nlev;ilev++)
    pops[ilev] /= sum;
}


/*....................................................................*/
void
getTransitionRates(struct transitionParams *user_data, double radius, double *Pops){
  
  int itemp,ipart,t_binlow,iline,k,l,ti, li, upper, lower, tnint=-1;
  double rnuc, Te, ne, aij, sigmaij, ve, bessel, ceij, gij, ceji;
  double interp_coeff, Qwater, vkin, collRate, tau, beta;
  
  molData *md = user_data->md;
  int ispec = user_data->ispec;
  configInfo *par = user_data->par;
  struct grid *gp = user_data->gp;
  int NEQ = user_data->array_size;
  double *A = user_data->A_array;
  
double *p = user_data->p_rates;
  double *coll_rates = user_data->coll_rates;

  /* 1. ONLY RECALCULATE COLLISIONS IF RADIUS HAS CHANGED (Optimization) */
  if (radius != user_data->last_radius) {
    double dens[md[ispec].npart], tkin[md[ispec].npart], LTEpops[md[ispec].nlev];
    
    rnuc = par->minScale;
    density(0.0,0.0,user_data->subGrid*radius,dens);
    temperature(0.0,0.0,user_data->subGrid*radius,tkin);
    lteOnePoint(md, ispec, tkin[0], LTEpops);
    
    double vel[3];
    velocity(0.,0.,user_data->subGrid*radius,vel);
    user_data->vexp = sqrt(vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2]);
    if (user_data->vexp < 1e-10) user_data->vexp = 1e-10; /* Prevent divide-by-zero */
    
    /* Initialize matrix with zeros */
    if(md[ispec].nlev<=0){
      if(!silent) bail_out("Matrix initialization error in solveStatEq");
      exit(1);
    }
    for(k=0; k < NEQ*NEQ; k++) coll_rates[k] = 0.0;
    
    /* Populate matrix with collisional transitions */
    for(ipart=0;ipart<md[ispec].npart;ipart++){
      struct cpData part = md[ispec].part[ipart];
      double *downrates = part.down;
      int di = md[ispec].part[ipart].densityIndex;
      if (di<0) continue;

      /* Collision temperature interpolation coefficients */
      if((tkin[ipart]>part.temp[0])&&(tkin[ipart]<part.temp[part.ntemp-1])){
        for(itemp=0;itemp<part.ntemp-1;itemp++){
          if((tkin[ipart]>part.temp[itemp])&&(tkin[ipart]<=part.temp[itemp+1])){
            tnint=itemp;
          }
        }
        interp_coeff =(tkin[ipart]-part.temp[tnint])/(part.temp[tnint+1]-part.temp[tnint]);
        t_binlow = tnint;
      } else if(tkin[ipart]<=part.temp[0]) {
        t_binlow = 0;
        interp_coeff = 0.0;
      } else {
        t_binlow = part.ntemp-2;
        interp_coeff = 1.0;
      }

      if (part.ntrans > 0){
        /* Use the LAMDA collision rates, if provided in the input file */ 
        for(ti=0;ti<part.ntrans;ti++){
          int coeff_index = ti*part.ntemp + t_binlow;
          double down = downrates[coeff_index] + interp_coeff*(downrates[coeff_index+1] - downrates[coeff_index]);
          double up = down*md[ispec].gstat[part.lcu[ti]]/md[ispec].gstat[part.lcl[ti]]
                    *exp(-HCKB*(md[ispec].eterm[part.lcu[ti]]-md[ispec].eterm[part.lcl[ti]])/tkin[ipart]);

          coll_rates[part.lcu[ti] * NEQ + part.lcl[ti]] += down*dens[ipart];
          coll_rates[part.lcl[ti] * NEQ + part.lcu[ti]] += up*dens[ipart];
        }
      } else {
        /* Otherwise use the Meudon approximation */
        vkin = sqrt(8.0*KBOLTZ*tkin[ipart]/PI * (1.0/md[ispec].amass + 1.0/MATM));
        collRate = dens[ipart]*vkin*XSEC*par->colliScale;
        /* CACHE OPTIMIZATION: Swapped l and k loops to prevent cache misses */
        for(l=0;l<md[ispec].nlev;l++){
          for(k=0;k<md[ispec].nlev;k++){
            coll_rates[l * NEQ + k] += collRate * LTEpops[k];  
          }
        }
      }
    }  
    
    /*GENERATE ELECTRON COLLISIONAL RATES (only for gas 0) AND ADD TO MATRIX*/
    //Presently, only electrons produced from collision parter 0 are considered.
    //It would be easy enough to add others, but the partner production rates would be needed as an input
    //parameter, like par->Qpartner, and their temperatures can be given as tkin[n] from temperature()
    /*Formalism of Zakharov et al. (2007)*/
    if(user_data->subGrid==SUBGRID1) Qwater = par->Q1;
    else Qwater = par->Q2;

    Te = Telec(radius,Qwater,tkin[0]); 
    ne = nelec(radius,Qwater,user_data->vexp,Te,par->rHelio,par->xne); 
     
    for(iline=0;iline<md[ispec].nline;iline++){
      aij = HPLANCK*md[ispec].freq[iline]/2./KBOLTZ/Te;
      sigmaij = ELEC_MASS*pow(ELEC_CHARGE,2)*pow(CLIGHT,3)*md[ispec].aeinst[iline]/16./pow(PI,2)/EPS_0/pow(HPLANCK,2)/pow(md[ispec].freq[iline],4);
      ve = sqrt(8.*KBOLTZ*Te/PI/ELEC_MASS);
      bessel = gsl_sf_bessel_K0(aij);
      ceij = ne*ve*sigmaij*2.*aij*exp(aij)*bessel;
      gij = md[ispec].gstat[md[ispec].lau[iline]]/md[ispec].gstat[md[ispec].lal[iline]];
      ceji = ne*ve*gij*sigmaij*2.*aij*exp(-aij)*bessel;

      coll_rates[md[ispec].lau[iline] * NEQ + md[ispec].lal[iline]] += ceij;
      coll_rates[md[ispec].lal[iline] * NEQ + md[ispec].lau[iline]] += ceji;
    }

    /* Add the pumping rates to the matrix */
    if(par->girdatfile!=NULL){
      /* CACHE OPTIMIZATION: Swapped l and k loops to prevent cache misses */
      for(l=0;l<md[ispec].nlev;l++){
        for(k=0;k<md[ispec].nlev;k++){
          if(k!=l) coll_rates[l * NEQ + k] += md[ispec].gir[l*md[ispec].nlev+k];
        }
      }
    }
    
    user_data->last_radius = radius;
  }

  /* 2. COPY CACHED RATES TO WORKING MATRIX */
  for(k=0; k < NEQ*NEQ; k++) p[k] = coll_rates[k];

  //Radiation trapping using the Escape Probability method
  if(par->useEP==1){ 
    double molDens[par->nSpecies];
    molNumDensity(0.0,0.0,user_data->subGrid*radius,molDens); 
    
    for(li=0;li<md[ispec].nline;li++){
      upper=md[ispec].lau[li];
      lower=md[ispec].lal[li];

      //Calculating the optical depth
      tau = ((A[li]*pow(CLIGHT,3))/(8*PI*pow(md[ispec].freq[li],3))) * 
            ((md[ispec].gstat[upper]/md[ispec].gstat[lower])*Pops[lower] - Pops[upper]) * 
            ((molDens[ispec]* radius)/user_data->vexp);

      //If the optical depth is small, ignore it
      if (tau > -1.0e-6 && tau < 1.0e-6){
        beta = 1.0;
      }else if (tau>0.0) {
        beta = (2/(3*tau)) - exp(-tau/2)*(tau*(gsl_sf_bessel_Kn(2,tau/2)-gsl_sf_bessel_K1(tau/2))/3 -  gsl_sf_bessel_K1(tau/2)); 
      }else if (tau<0.0){ 
        beta = (1 - exp(-tau)) / tau;
      }
      p[upper * NEQ + lower] += A[li]*beta;
    }//end for
  }//end if
  //No photon trapping simulation
  else if(par->useEP==0){
    for(li=0;li<md[ispec].nline;li++){
      upper=md[ispec].lau[li];
      lower=md[ispec].lal[li];
      p[upper * NEQ + lower] += A[li];
    }
  }
  //Calculate photon trapping self-consistently, using full 3D treatment of radiation field, and photon propagation
  else if(par->useEP==2){
    int gp_pIntensity = user_data->gp_pIntensity;
    int *gp_sorter = user_data->gp_sorter;
    double *jbar_grid = user_data->jbar_grid;
    
    if(radius >= gp[gp_sorter[gp_pIntensity-1]].radius){
      for(li=0;li<md[ispec].nline;li++){
        upper=md[ispec].lau[li];
        lower=md[ispec].lal[li];
        // Note - this is broken as jbar_grid should refer to the current subgrid, not the total grid, and we should be choosing the closest grid point to the cone center, relative to our present radius     
        p[upper * NEQ + lower] += md[ispec].beinstl[li]*jbar_grid[(par->pIntensity-1) * md[ispec].nline + li];
        p[lower * NEQ + upper] += md[ispec].beinstl[li]*jbar_grid[(par->pIntensity-1) * md[ispec].nline + li];
      }
    } else {
      for(li=0;li<md[ispec].nline;li++){
        upper=md[ispec].lau[li];
        lower=md[ispec].lal[li];
        // Note - this calculation is badly broken as jbar_grid should refer to the current subgrid, not the total grid, and we should be choosing the closest grid point to the cone center, relative to our present radius... (par->pIntensity-1) is used only as a placeholder grid point, to allow the code to run without causing major problems
        p[upper * NEQ + lower] += md[ispec].beinstl[li]*jbar_grid[(par->pIntensity-1) * md[ispec].nline + li];
        p[lower * NEQ + upper] += md[ispec].beinstl[li]*jbar_grid[(par->pIntensity-1) * md[ispec].nline + li];
      }
    }
  }
}

/*....................................................................*/
void
LTE(configInfo *par, struct grid *gp, molData *md){
  int id,ispec;

  for(id=0;id<par->pIntensity;id++){
    for(ispec=0;ispec<par->nSpecies;ispec++){
      lteOnePoint(md, ispec, gp[id].t[0], gp[id].mol[ispec].pops);
    }
  }

}


/*....................................................................*/
/* Sets the Differential Equation (Pdot) to be solved by CVode */
int f(realtype radius, N_Vector P, N_Vector Pdot, void *data){

  int i, j, NEQ;
  struct transitionParams *user_data = data;
  NEQ = user_data->array_size;
  double *p = user_data->p_rates;
  double Pops_array[NEQ];

  // Passing P values to Pops_array for readability
  for(i=0; i < NEQ; ++i)
    Pops_array[i] = Ith(P,i);

  // Call the updated function
  getTransitionRates(user_data, radius, Pops_array);

  // Initialize Pdot
  for(i=0; i < NEQ; ++i)
    Ith(Pdot,i) = 0.0;

  // Cache-friendly matrix-vector product
  
  // 1. Add population coming FROM j TO i
  for(j=0; j < NEQ; ++j){
    double pop_j = Pops_array[j];
    for(i=0; i < NEQ; ++i){
      if(i != j) {
        Ith(Pdot,i) += pop_j * p[j * NEQ + i];
      }
    }
  }

  // 2. Subtract population going FROM i TO j
  for(i=0; i < NEQ; ++i){
    double pop_i = Pops_array[i];
    double out_rate_sum = 0.0;
    for(j=0; j < NEQ; ++j){
      if(i != j) {
        out_rate_sum += p[i * NEQ + j];
      }
    }
    Ith(Pdot,i) -= pop_i * out_rate_sum;
  }
   
  // Multiply equations by 1/V to get dP/dr     
  double inv_vexp = 1.0 / user_data->vexp;
  for(i=0; i < NEQ; ++i){
     Ith(Pdot,i) *= inv_vexp;
  }
 
  return(0);
}

/*....................................................................*/
/* reduceTol: NOW TAKES reltol BY REFERENCE.
   Previously reltol was passed by value, so the *reltol in the caller
   (solveStatEq) never actually decreased - only abstol (passed by pointer)
   did. This meant the surrounding while(cvstatus<0 && reltol>MINTOL) loop
   in solveStatEq checked a stale value and risked looping indefinitely.
   The caller now does "reduceTol(abstol, &reltol, ...)" and no longer
   needs its own separate "reltol *= 0.1" afterwards. */
void reduceTol(N_Vector abstol, realtype *reltol, void *cvode_mem, double factor, int NEQ){
   printf("INFO: Reducing RTOL and ATOL by 0.1 and restarting CVODE\n");
   int i, retval;

   for(i=0; i < NEQ; ++i){
      Ith(abstol,i) = Ith(abstol,i) * factor;
   }

   *reltol = (*reltol) * factor;

   retval = CVodeSVtolerances(cvode_mem, *reltol, abstol);
   if (check_retval(&retval, "CVodeSVtolerances", 1)) return;
}

/*....................................................................*/
/* Exact Jacobian-times-vector product, for use with CVodeSetJacTimes.

   VALID when par->useEP == 0 OR par->useEP == 2: in both cases the rate
   matrix p does not depend on the current population vector Pops (for
   useEP==2, p depends only on jbar_grid and radius, not Pops), so
   f(t,P) is exactly linear: f = (1/vexp)*p^T*P, and J = (1/vexp)*p^T
 */

int JacTimesVec(N_Vector v, N_Vector Jv, realtype t, N_Vector P, N_Vector fP,
                 void *user_data_ptr, N_Vector tmp){

  struct transitionParams *user_data = user_data_ptr;
  int NEQ   = user_data->array_size;
  double *p = user_data->p_rates;
  int i,j;

  for(i=0; i < NEQ; ++i)
    Ith(Jv,i) = 0.0;

  /* Contribution flowing FROM j TO i */
  for(j=0; j < NEQ; ++j){
    double v_j = Ith(v,j);
    for(i=0; i < NEQ; ++i){
      if(i != j)
        Ith(Jv,i) += v_j * p[j * NEQ + i];
    }
  }

  /* Contribution flowing FROM i TO j (subtracted) */
  for(i=0; i < NEQ; ++i){
    double out_rate_sum = 0.0;
    for(j=0; j < NEQ; ++j){
      if(i != j)
        out_rate_sum += p[i * NEQ + j];
    }
    Ith(Jv,i) -= Ith(v,i) * out_rate_sum;
  }

  double inv_vexp = 1.0 / user_data->vexp;
  for(i=0; i < NEQ; ++i)
    Ith(Jv,i) *= inv_vexp;

  return 0;
}

/*....................................................................*/
void
solveStatEq(struct grid *gp, molData *md, const int ispec, configInfo *par\
  , struct blendInfo blends, int *nextMolWithBlend, gridPointData **mp\
  , double **halfFirstDs, int *nMaserWarnings,struct grid *gp3D,int gp_pIntensity, int gp_ncell, double *radii, int subGrid_pIntensity, int *gp_sorter, int subGrid){
  
  int id;
  realtype reltol, t;
  N_Vector P, abstol;
  void *cvode_mem;
  int i, j, k, cvodeErrs = 0;
  int retval, cvstatus, index;
  
  /* THREAD SAFETY FIX: Create a local copy of par to prevent Race Conditions */
  configInfo local_par = *par;
  
  int ncell = local_par.ncell, pIntensity = local_par.pIntensity;

  /* Initializing parameters to be used by CVode */
  int NEQ = md[ispec].nlev; // Number of Equations
  double A[md[ispec].nline]; // Einstein As
  double Pops[NEQ]; // Level Populations
  double (*popGrid)[NEQ]; // Populations as a function of radius 
  double (*jbar_grid)[md[ispec].nline] = NULL;
  
  popGrid = malloc(sizeof(double[NRADS][NEQ]));
  
  for(id=0;id<md[ispec].nline;id++)
    A[id] = md[ispec].aeinst[id];
  
  if (local_par.useEP==2){
    jbar_grid = malloc(sizeof(double[pIntensity][md[ispec].nline])); 
  
    // Changing the values to those of the full grid to perform the jbar update
    local_par.ncell = gp_ncell;
    local_par.pIntensity = gp_pIntensity; 

    for(id=0;id<local_par.pIntensity;id++){
      updateJBar(id, md, gp3D, ispec, &local_par, blends, nextMolWithBlend[id], mp[id], halfFirstDs[id]);
    }
    for(i=0;i<local_par.pIntensity;i++){
        // OPTIMIZATION: Hoist k_match search out of the j-loop
        int k_match = 0;
        for(k=0;k<gp_pIntensity;k++){
            if(gp3D[i].id==gp[k].id){ 
                k_match = k;
                break;
            }
        }
        for(j=0;j<md[ispec].nline;j++){
            jbar_grid[k_match][j] = mp[i][ispec].jbar[j];
        }
    }
    // Reverting to the original values
    local_par.ncell = ncell;
    local_par.pIntensity = pIntensity;
  }

  /* Initializing Pops */
  for(i=0;i<md[ispec].nlev;i++){
    Pops[i] = gp[gp_sorter[0]].mol[ispec].pops[i]; 
    popGrid[0][i] = Pops[i];
  }

  /* ALLOCATE CACHED MATRICES ONCE */
  double *p_rates = malloc(sizeof(double) * NEQ * NEQ);
  double *coll_rates = malloc(sizeof(double) * NEQ * NEQ);

  /* Pass the local_par into user_data to maintain thread safety */
  struct transitionParams user_data = {
      NEQ, A, md, ispec, gp, &local_par, (double *)jbar_grid, nMaserWarnings, gp_sorter, 
      subGrid, subGrid_pIntensity, 
      p_rates, coll_rates, -1.0, 1.0 /* initializers for caching */
  }; 

  P = abstol = NULL;
  cvode_mem = NULL;
  
  /* Fallback solver variables */
  SUNLinearSolver LS = NULL;
  SUNMatrix sunMatrix = NULL;
  int using_dense; 

  /* Create serial vector of length NEQ for I.C. and abstol */
  P = N_VNew_Serial(NEQ);
  if (check_retval((void *)P, "N_VNew_Serial", 0)) return;
  abstol = N_VNew_Serial(NEQ); 
  if (check_retval((void *)abstol, "N_VNew_Serial", 0)) return;

  /* Set the scalar relative tolerance */
  reltol = RTOL;

  for(i=0; i < NEQ; ++i){
    Ith(P,i) = Pops[i];
    Ith(abstol,i) = ATOL;
  }

  cvode_mem = CVodeCreate(CV_BDF);
  if (check_retval((void *)cvode_mem, "CVodeCreate", 0)) return;

  retval = CVodeInit(cvode_mem, f, radii[0], P);
  if (check_retval(&retval, "CVodeInit", 1)) return;

  retval = CVodeSetUserData(cvode_mem, &user_data);
  if (check_retval(&retval, "CVodeSetUserData", 0)) return;

  retval = CVodeSVtolerances(cvode_mem, reltol, abstol);
  if (check_retval(&retval, "CVodeSVtolerances", 1)) return;

  retval = CVodeSetMaxNumSteps(cvode_mem, 5000);
  if (check_retval(&retval, "CVodeSetMaxNumSteps", 1)) return;

  /* --- SOLVER SELECTION ---
     Force Dense whenever useEP==1: the escape-probability term makes p
     depend on Pops (nonlinear system), invalidating JacTimesVec. Note this
     is deliberately useEP==1 specifically, NOT useEP>0 - useEP==2 (full 3D
     photon trapping via jbar_grid) does NOT depend on Pops and is safe for
     SPGMR+exact-Jacobian, same as useEP==0. For useEP==0 or useEP==2,
     prefer SPGMR above DENSE_NEQ_CUTOFF, backed by the exact
     Jacobian-vector product and a widened Krylov subspace (20, rather
     than SUNDIALS' small default) to guard against the convergence
     failures that made SPGMR's earlier (FD-Jacobian) behaviour unreliable. */
  int forceDense = (local_par.useEP == 1);

  if (!forceDense && NEQ > DENSE_NEQ_CUTOFF) {
    using_dense = 0;
    LS = SUNLinSol_SPGMR(P, PREC_NONE, 20);
    if(check_retval((void *)LS, "SUNLinSol_SPGMR", 0)) return;

    retval = CVodeSetLinearSolver(cvode_mem, LS, NULL);
    if(check_retval(&retval, "CVodeSetLinearSolver", 1)) return;

    retval = CVodeSetJacTimes(cvode_mem, NULL, JacTimesVec);
    if(check_retval(&retval, "CVodeSetJacTimes", 1)) return;
  } else {
    using_dense = 1;
    sunMatrix = SUNDenseMatrix(NEQ, NEQ);
    if(check_retval((void *)sunMatrix, "SUNDenseMatrix", 0)) return;

    LS = SUNLinSol_Dense(P, sunMatrix);
    if(check_retval((void *)LS, "SUNLinSol_Dense", 0)) return;

    retval = CVodeSetLinearSolver(cvode_mem, LS, sunMatrix);
    if(check_retval(&retval, "CVodeSetLinearSolver", 1)) return;
  }

  printf("Linear solver selected for sub-grid %d, species %d: %s (NEQ=%d, useEP=%d)\n",
         subGrid, ispec, using_dense ? "Dense" : "Krylov", NEQ, local_par.useEP);
  fflush(stdout);

  retval = CVodeSetJacFn(cvode_mem, NULL);
  if(check_retval(&retval, "CVodeSetJacFn", 1)) return;
  
  printf("Starting CVODE time-dependent solver for sub-grid %d, species %d...\n",subGrid,ispec);
  fflush(stdout);

  // Call CVODE for each radius
  do{
     for(i=1; i<NRADS; i++){
         cvstatus = CVode(cvode_mem, radii[i], P, &t, CV_NORMAL);
        
       if(cvstatus == CV_SUCCESS){
         for(j=0;j<md[ispec].nlev;j++){ 
           popGrid[i][j] = Ith(P,j); 
         }
       }
      
       /* IF SOLVER FAILS */
       if(cvstatus < 0 && reltol > MINTOL) { 
         
         /* Reset Populations to initial values before restarting */
         for(k=0; k < NEQ; ++k) Ith(P,k) = popGrid[0][k];
         
         /* TIER 1 FALLBACK: SWITCH TO DENSE SOLVER */
         if (using_dense == 0) {
             printf("WARNING: SPGMR failed (stiff system). Switching to Dense Direct Solver...\n");
             using_dense = 1;
             
             if (LS != NULL) SUNLinSolFree(LS);
             
             retval = CVodeInit(cvode_mem, f, radii[0], P);
             CVodeSetUserData(cvode_mem, &user_data);
             CVodeSetMaxNumSteps(cvode_mem, 5000);
             
             sunMatrix = SUNDenseMatrix(NEQ, NEQ);
             LS = SUNLinSol_Dense(P, sunMatrix);
             CVodeSetLinearSolver(cvode_mem, LS, sunMatrix);
             CVodeSetJacFn(cvode_mem, NULL); 
             
             break; /* Break spatial loop to restart integration */
         } 
         /* TIER 2 FALLBACK: DENSE SOLVER FAILED, REDUCE TOLERANCES */
         else {
             printf("WARNING: Dense solver failed. Reducing tolerances...\n");
             
             retval = CVodeInit(cvode_mem, f, radii[0], P);
             CVodeSetUserData(cvode_mem, &user_data);
             CVodeSetMaxNumSteps(cvode_mem, 5000);
             
             if (LS != NULL) SUNLinSolFree(LS);
             if (sunMatrix != NULL) SUNMatDestroy(sunMatrix);
             
             sunMatrix = SUNDenseMatrix(NEQ, NEQ);
             LS = SUNLinSol_Dense(P, sunMatrix);
             CVodeSetLinearSolver(cvode_mem, LS, sunMatrix);
             CVodeSetJacFn(cvode_mem, NULL);
             
             /* reduceTol now takes reltol by reference and updates it
                correctly in place - no separate "reltol *= 0.1" needed. */
             reduceTol(abstol, &reltol, cvode_mem, 0.1, NEQ);
             
             break; /* Break spatial loop to restart integration */
         }
      } 
      else if (cvstatus > 0) {
          cvodeErrs++;
          printf("CVODE error %d (continuing to next timestep - check populations!!)\n",cvstatus);
      }
      
      if(cvodeErrs >= 15){ 
         bail_out("CVODE solver failure - check physical model and tolerances.");
         exit(1);
      }
    }
  } while (cvstatus < 0 && reltol > MINTOL);

  fflush(stdout);
  
  /* Interpolate the computed radial populations onto the Delaunay grid */ 
  for(i=0;i<subGrid_pIntensity;i++){
     index = NRADS - 1; /* Default to outer boundary to prevent memory read violations */
     for(j=1;j<NRADS;j++){ 
        if(radii[j]>gp[i].radius){
           index = j;
           break;
        }
     }   
     for(k=0;k<md[ispec].nlev;k++){
        gp[i].mol[ispec].pops[k] = linear_interp(radii[index-1],radii[index],popGrid[index-1][k],popGrid[index][k],gp[i].radius); 
     }
  }

  /* Free CVODE Vectors and Memory */
  N_VDestroy(P);
  N_VDestroy(abstol);
  CVodeFree(&cvode_mem);

  /* Free the linear solver & matrix memory */
  if (LS != NULL) SUNLinSolFree(LS);
  if (sunMatrix != NULL) SUNMatDestroy(sunMatrix);
  
  /* Free cached arrays */
  free(p_rates);
  free(coll_rates);
  free(popGrid);
  
  if (local_par.useEP==2)
    free(jbar_grid);
}

/*....................................................................*/
int
levelPops(molData *md, configInfo *par, struct grid *gp, int *popsdone, double *lamtab, double *kaptab, const int nEntries, struct grid *gp3D, int gp_pIntensity, int gp_ncell, double *radii, int subGrid_pIntensity, int *gp_sorter, int subGrid){

  int id,ispec,i,nVerticesDone,nlinetot;
  int ncell = par->ncell, pIntensity = par->pIntensity;
  int totalNMaserWarnings=0;
  const gsl_rng_type *ranNumGenType = gsl_rng_ranlxs2;
  struct blendInfo blends;
  char message[STR_LEN_0];
  int RNG_seeds[par->nThreads];
  gsl_error_handler_t *defaultErrorHandler=NULL;
  int nextMolWithBlend[gp_pIntensity],nMaserWarnings[gp_pIntensity];
  for(id=0;id<gp_pIntensity;id++){
    nMaserWarnings[id] = 0;
  }

  nlinetot = 0;
  for(ispec=0;ispec<par->nSpecies;ispec++){
    nlinetot += md[ispec].nline;
  }  

  if(par->lte_only){
    LTE(par,gp,md);
    if(par->outputfile) popsout(par,gp,md);
   
  }else{
  /* Non-LTE */
  
   /* Random number generator */
    gsl_rng *ran = gsl_rng_alloc(ranNumGenType);
    if(fixRandomSeeds)
      gsl_rng_set(ran, 1237106) ;
    else 
      gsl_rng_set(ran,time(0));

    gsl_rng **threadRans;
    threadRans = malloc(sizeof(gsl_rng *)*gp_pIntensity);

    for (i=0;i<gp_pIntensity;i++){
      threadRans[i] = gsl_rng_alloc(ranNumGenType);
      if (par->resetRNG==1) RNG_seeds[i] = (int)(gsl_rng_uniform(ran)*1e6);
      else gsl_rng_set(threadRans[i],(int)(gsl_rng_uniform(ran)*1e6));
    }

   calcGridCollRates(par,md,gp);
   freeGridCont(par, gp);
   mallocGridCont(par, md, gp);
   calcGridLinesDustOpacity(par, md, lamtab, kaptab, nEntries, gp);

   /* Check for blended lines */
   lineBlend(md, par, &blends);

   /* Initialize populations with Boltzmann distribution (assuming LTE) */
   if (par->useEP != 2) LTE(par,gp,md); //If useEp ==2, then the populations have already been initialized, and we don't want to override them

   defaultErrorHandler = gsl_set_error_handler_off();

    gridPointData *mp[gp_pIntensity];
    double *halfFirstDs[gp_pIntensity];

    calcGridMolSpecNumDens(par,md,gp);
    totalNMaserWarnings = 0;
    nVerticesDone=0;

    //TODO: This for loop could be parallelized
    for(id=0;id<gp_pIntensity;id++){
      ++nVerticesDone;
      nMaserWarnings[id]=0;
      nextMolWithBlend[id] = 0;
      mp[id]=malloc(sizeof(gridPointData)*par->nSpecies);
      halfFirstDs[id] = malloc(sizeof(*halfFirstDs)*gp3D[id].nphot);

      for (ispec=0;ispec<par->nSpecies;ispec++){
        mp[id][ispec].jbar = malloc(sizeof(double)*md[ispec].nline);
        mp[id][ispec].phot = malloc(sizeof(double)*md[ispec].nline*gp3D[id].nphot);
        mp[id][ispec].vfac = malloc(sizeof(double)*                gp3D[id].nphot);
        mp[id][ispec].vfac_loc = malloc(sizeof(double)*            gp3D[id].nphot);

      }
      if(gp3D[id].dens[0] < 0 && gp3D[id].t[0] < 0){
        printf("\nError on grid point = %d\n, aborting", id);
        exit(1);
      }
      else if (par->useEP==2){
        par->pIntensity = gp_pIntensity;
        par->ncell  = gp_ncell;
        calculateJBar(id,gp3D,md,threadRans[id],par,nlinetot,blends,mp[id],halfFirstDs[id],&nMaserWarnings[id]);
      }
    }

    //These corresponds to the values of the subgrid, which would have been previously changed to the values that correspond to the complete grid if useEP=2 for the jbar calculation
    par->pIntensity = pIntensity;
    par->ncell = ncell;

    for(ispec=0;ispec<par->nSpecies;ispec++){
      solveStatEq(gp,md,ispec,par,blends,nextMolWithBlend,mp,halfFirstDs, nMaserWarnings,gp3D,gp_pIntensity, gp_ncell, radii, subGrid_pIntensity, gp_sorter, subGrid); 
      for(i=0;i<gp_pIntensity;i++)
        if(par->blend && blends.mols!=NULL && ispec==blends.mols[nextMolWithBlend[i]].molI)
          nextMolWithBlend[i] = nextMolWithBlend[i] + 1;
    }

    printf("SolveStatEq: DONE\n");
    fflush(stdout);

    for(id=0;id<gp_pIntensity;id++){
      totalNMaserWarnings += nMaserWarnings[id];
    }

    if(!silent && totalNMaserWarnings>0){
      snprintf(message, STR_LEN_0, "Maser warning: optical depth dropped below -%4.1f %d times this iteration.", MAX_NEG_OPT_DEPTH, totalNMaserWarnings);
      warning(message);
    }

    if(!silent) warning("");
    for (i=0;i<gp_pIntensity;i++){
      freeGridPointData(par->nSpecies, mp[i]);
      free(halfFirstDs[i]);
    }

  freeMolsWithBlends(blends.mols, blends.numMolsWithBlends);
  for (i=0;i<par->pIntensity;i++)
    gsl_rng_free(threadRans[i]);
  free(threadRans);
  gsl_rng_free(ran);
 }//end else

  par->dataFlags |= (1 << DS_bit_populations);

  if(par->binoutputfile != NULL) binpopsout(par,gp,md);

  *popsdone=1;

  return (1);
}