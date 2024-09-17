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
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
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
  double *transition_rates; //Holds transition rates
  molData *md;
  int ispec;
  struct grid *gp;
  configInfo *par;
  double *jbar_grid;
  int *nMaserWarnings;
  int *gp_sorter;
  int subGrid;
  int gp_pIntensity;
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
  }
}

/*....................................................................*/
void lineBlend(molData *m, configInfo *par, struct blendInfo *blends){
  /*
This obtains information on all the lines of all the radiating species which have other lines within some cutoff velocity separation.

A variable of type 'struct blendInfo' has a nested structure which can be illustrated diagrammaticaly as follows.

  Structs:  blendInfo   molWithBlends   lineWithBlends    blend

  Variables:  blends
      .numMolsWithBlends     ____________________
      .*mols--------------->|.molI               |
                            |.numLinesWithBlends |   ___________
                            |.*lines--------------->|.lineI     |
                            |____________________|  |.numBlends |           ________
                            |        etc         |  |.*blends------------->|.molJ   |
                                                    |___________|          |.lineJ  |
                                                    |    etc    |          |.deltaV |
                                                                           |________|
                                                                           |   etc  |

Pointers are indicated by a * before the attribute name and an arrow to the memory location pointed to.
  */
  int molI, lineI, molJ, lineJ;
  int nmwb, nlwb, numBlendsFound, li, bi;
  double deltaV;
  struct blend *tempBlends=NULL;
  struct lineWithBlends *tempLines=NULL;

  /* Dimension blends.mols first to the total number of species, then realloc later if need be.
  */
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
The idea here is to select for the next grid point, that one which lies closest (with a little randomizing jitter) to the photon track, while requiring the direction of the edge to be in the 'forward' hemisphere of the photon direction.

Note that this is called from within the multi-threaded block.
  */
  int i,ni,niOfSmallest=-1,niOfNextSmallest=-1;
  double dirCos,distAlongTrack,dirFromStart[3],coord,distToTrackSquared,smallest=0.0,nextSmallest=0.0;
  const static double scatterReduction = 0.4;
  /*
This affects the ratio of N_2/N_1, where N_2 is the number of times the edge giving the 2nd-smallest distance from the photon track is chosen and N_1 ditto the smallest. Some ratio values obtained from various values of scatterReduction:

  scatterReduction  <N_2/N_1>
    1.0     0.42
    0.5     0.75
    0.4     0.90
    0.2     1.52

Note that the equivalent ratio value produced by the 1.6 code was 0.91.
  */

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

  /* Choose the edge to follow.
  */
  if(i>1){ /* then nextSmallest, niOfNextSmallest should exist. */
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
  /*
Note that this is called from within the multi-threaded block.
  */

  /* convolution of a Gaussian with a box */
  double binv_this, binv_next, v[5];

  binv_this=g[id].mol[molI].binv;
  binv_next=(g[id].neigh[k])->mol[molI].binv;
  v[0]=deltav-dotProduct3D(inidir,g[id].vel);
  v[1]=deltav-dotProduct3D(inidir,&(g[id].v1[3*k]));
  v[2]=deltav-dotProduct3D(inidir,&(g[id].v2[3*k]));
  v[3]=deltav-dotProduct3D(inidir,&(g[id].v3[3*k]));
  v[4]=deltav-dotProduct3D(inidir,g[id].neigh[k]->vel);

  /* multiplying by the appropriate binv changes from velocity to doppler widths(?) */
  /* if the values were be no more than 2 erf table bins apart, we just take a single Gaussian */

  /*
  vfac_out is the lineshape for the part of the edge in the current Voronoi cell,
  vfac_in is for the part in the next cell
  */

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
  /*
Note that this is called from within the multi-threaded block.
  */

  /* convolution of a Gaussian with a box */
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

    /* Choose random initial photon direction (the distribution used here is even over the surface of a sphere of radius 1).
    */
    pt_theta=gsl_rng_uniform(ran)*2*M_PI;
    pt_z=2*gsl_rng_uniform(ran)-1;
    semiradius = sqrt(1.-pt_z*pt_z);
    inidir[0]=semiradius*cos(pt_theta);
    inidir[1]=semiradius*sin(pt_theta);
    inidir[2]=pt_z;

    /* Choose the photon frequency/velocity offset.
    */
    segment=gsl_rng_uniform(ran)-0.5;
    /*
    Values of segment should be evenly distributed (considering the
    entire ensemble of photons) between -0.5 and +0.5.
    */

    for (molI=0;molI<par->nSpecies;molI++){
      /* Is factor 4.3=[-2.15,2.15] enough?? */
      deltav[molI]=4.3*segment*gp[id].mol[molI].dopb+dotProduct3D(inidir,gp[id].vel);
      /*
      This is the local (=evaluated at a grid point, not averaged over the local cell) lineshape.
      We store this for later use in ALI loops.
      */
      mp[molI].vfac_loc[iphot]=gaussline(deltav[molI]-dotProduct3D(inidir,gp[id].vel),gp[id].mol[molI].binv);
    }

    here = gp[id].id;

    /* Photon propagation loop */
    numLinks=0;
    while(!gp[here].sink){ /* Testing for sink at loop start is redundant for the first step, since we only start photons from non-sink points, but it makes for simpler code. */
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
        /*
        Contribution of the local cell to emission and absorption is done in updateJBar.
        We only store the vfac for the local cell for use in ALI loops.
        */
        here=there;
    continue;
      }

      /* If we've got to here, we have progressed beyond the first edge. Length of the new "in" edge is the length of the previous "out".
      */
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

          /* cont and blend could use the same alpha and jnu counter, but maybe it's clearer this way */

          /* Line blending part.
          */
          if(par->blend && blends.mols!=NULL && molI==blends.mols[nextMolWithBlend].molI\
          && lineI==blends.mols[nextMolWithBlend].lines[nextLineWithBlend].lineI){

            for(bi=0;bi<blends.mols[nextMolWithBlend].lines[nextLineWithBlend].numBlends;bi++){
              molJ  = blends.mols[nextMolWithBlend].lines[nextLineWithBlend].blends[bi].molJ;
              lineJ = blends.mols[nextMolWithBlend].lines[nextLineWithBlend].blends[bi].lineJ;
              velProj = deltav[molI] - blends.mols[nextMolWithBlend].lines[nextLineWithBlend].blends[bi].deltaV;
        /*  */
              if(par->edgeVelsAvailable)
                calcLineAmpPWLin(gp,here,neighI,molJ,velProj,inidir,&vblend_in,&vblend_out);
              else
                calcLineAmpLin(gp,here,neighI,molJ,velProj,inidir,&vblend_in,&vblend_out);

        /* we should use also the previous vblend_in, but I don't feel like writing the necessary code now */
              sourceFunc_line(&md[molJ],vblend_out,&(gp[here].mol[molJ]),lineJ,&jnu_blend,&alpha_blend);
              /* note that sourceFunc* increment jnu and alpha, they don't overwrite it  */
            }

            nextLineWithBlend++;
            if(nextLineWithBlend>=blends.mols[nextMolWithBlend].numLinesWithBlends){
              nextLineWithBlend = 0;
              /* The reason for doing this is as follows. Firstly, we only enter the present IF block if molI has at least 1 line which is blended with others; and further, if we have now processed all blended lines for that molecule. Thus no matter what value lineI takes for the present molecule, it won't appear as blends.mols[nextMolWithBlend].lines[i].lineI for any i. Yet we will still test blends.mols[nextMolWithBlend].lines[nextLineWithBlend], thus we want nextLineWithBlend to at least have a sensible value between 0 and blends.mols[nextMolWithBlend].numLinesWithBlends-1. We could set nextLineWithBlend to any number in this range in safety, but zero is simplest. */
            }
          }
          /* End of line blending part */

    /* as said above, out-in split should be done also for blended lines... */

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
        } /* Next line this molecule. */

        if(par->blend && blends.mols!=NULL && molI==blends.mols[nextMolWithBlend].molI)
          nextMolWithBlend++;
      }

      here=there;
    };

    /* Add cmb contribution.
    */
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

        /* Line blending part.
        */
        if(par->blend && blends.mols!=NULL && molI==blends.mols[nextMolWithBlend].molI\
        && lineI==blends.mols[nextMolWithBlend].lines[nextLineWithBlend].lineI){
          for(bi=0;bi<blends.mols[nextMolWithBlend].lines[nextLineWithBlend].numBlends;bi++){
            molJ  = blends.mols[nextMolWithBlend].lines[nextLineWithBlend].blends[bi].molJ;
            lineJ = blends.mols[nextMolWithBlend].lines[nextLineWithBlend].blends[bi].lineJ;
            /*
            The next line is not quite correct, because vfac may be different for other molecules due to different values of binv. Unfortunately we don't necessarily have vfac for molJ available yet.
            */
            sourceFunc_line(&md[molJ],mp[molI].vfac[iphot],&(gp[posn].mol[molJ]),lineJ,&jnu,&alpha);
      /* note that sourceFunc* increment jnu and alpha, they don't overwrite it  */
          }

          nextLineWithBlend++;
          if(nextLineWithBlend>=blends.mols[nextMolWithBlend].numLinesWithBlends){
            nextLineWithBlend = 0;
            /* The reason for doing this is as follows. Firstly, we only enter the present IF block if molI has at least 1 line which is blended with others; and further, if we have now processed all blended lines for that molecule. Thus no matter what value lineI takes for the present molecule, it won't appear as blends.mols[nextMolWithBlend].lines[i].lineI for any i. Yet we will still test blends.mols[nextMolWithBlend].lines[nextLineWithBlend], thus we want nextLineWithBlend to at least have a sensible value between 0 and blends.mols[nextMolWithBlend].numLinesWithBlends-1. We could set nextLineWithBlend to any number in this range in safety, but zero is simplest. */
          }
        }
        /* End of line blending part */

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
getTransitionRates(molData *md, int ispec, struct grid *gp, configInfo *par, int NEQ, double A[NEQ-1], double *p, double radius, double *jbar_grid, double *Pops, int *nMaserWarnings, int *gp_sorter, int subGrid, int subGrid_pIntensity, double vexp){
  int itemp,ipart,t_binlow,iline,k,l,ti, li, upper, lower, j,tnint=-1;
  double rnuc, Te, ne, aij, sigmaij, ve, bessel, ceij, gij, ceji, dens[md[ispec].npart], tkin[md[ispec].npart], LTEpops[md[ispec].nlev];
  double jbar[par->pIntensity], interp_coeff, Qwater;
  double vkin,collRate;

  rnuc = par->minScale;
  density(0.0,0.0,subGrid*radius,dens);
  temperature(0.0,0.0,subGrid*radius,tkin);
  lteOnePoint(md, ispec, tkin[0], LTEpops);

  /* Initialize matrix with zeros */
  if(md[ispec].nlev<=0){
    if(!silent) bail_out("Matrix initialization error in solveStatEq");
    exit(1);
  }

  for(k=0;k<md[ispec].nlev;k++)
    for(l=0;l<md[ispec].nlev;l++)
      p[k * NEQ + l] = 0.0;

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
      double down = downrates[coeff_index]\
                  + interp_coeff*(downrates[coeff_index+1] - downrates[coeff_index]);
      double up = down*md[ispec].gstat[part.lcu[ti]]/md[ispec].gstat[part.lcl[ti]]\
                *exp(-HCKB*(md[ispec].eterm[part.lcu[ti]]-md[ispec].eterm[part.lcl[ti]])/tkin[ipart]);

      p[part.lcu[ti] * NEQ + part.lcl[ti]] = p[part.lcu[ti] * NEQ + part.lcl[ti]] + down*dens[ipart];
      p[part.lcl[ti] * NEQ + part.lcu[ti]] = p[part.lcl[ti] * NEQ + part.lcu[ti]] + up*dens[ipart];
    }
   }else{
   /* Otherwise use the Meudon approximation */
     vkin = sqrt(8.0*KBOLTZ*tkin[ipart]/PI * (1.0/md[ispec].amass + 1.0/MATM));
     collRate = dens[ipart]*vkin*XSEC*par->colliScale;
     for(k=0;k<md[ispec].nlev;k++){
        for(l=0;l<md[ispec].nlev;l++){
          p[l * NEQ + k] = p[l * NEQ + k] + collRate * LTEpops[k];  
        }
     }

   }

  }  
  /*GENERATE ELECTRON COLLISIONAL RATES (only for gas 0) AND ADD TO MATRIX*/
  //Presently, only electrons produced from collision parter 0 are considered.
  //It would be easy enough to add others, but the partner production rates would be needed as an input
  //parameter, like par->Qpartner, and their temperatures can be given as tkin[n] from temperature()
  /*Formalism of Zakharov et al. (2007)*/

   if(subGrid==SUBGRID1){
    Qwater = par->Q1;
   }else{
    Qwater = par->Q2;
   }

   Te = Telec(radius,Qwater,tkin[0]); 
   ne = nelec(radius,Qwater,vexp,Te,par->rHelio,par->xne); 
   
   for(iline=0;iline<md[ispec].nline;iline++){
     aij = HPLANCK*md[ispec].freq[iline]/2./KBOLTZ/Te;
     sigmaij = ELEC_MASS*pow(ELEC_CHARGE,2)*pow(CLIGHT,3)*md[ispec].aeinst[iline]/16./pow(PI,2)/EPS_0/pow(HPLANCK,2)/pow(md[ispec].freq[iline],4);
     ve = sqrt(8.*KBOLTZ*Te/PI/ELEC_MASS);
     bessel = gsl_sf_bessel_K0(aij);
     ceij = ne*ve*sigmaij*2.*aij*exp(aij)*bessel;
     gij = md[ispec].gstat[md[ispec].lau[iline]]/md[ispec].gstat[md[ispec].lal[iline]];
     ceji = ne*ve*gij*sigmaij*2.*aij*exp(-aij)*bessel;

     p[md[ispec].lau[iline] * NEQ + md[ispec].lal[iline]] = p[md[ispec].lau[iline] * NEQ + md[ispec].lal[iline]] + ceij;
     p[md[ispec].lal[iline] * NEQ + md[ispec].lau[iline]] = p[md[ispec].lal[iline] * NEQ + md[ispec].lau[iline]] + ceji;

   }

   /* Add the pumping rates to the matrix */
   for(k=0;k<md[ispec].nlev;k++)
    for(l=0;l<md[ispec].nlev;l++)
      if(par->girdatfile!=NULL)
        if(k!=l){
          if(par->girdatfile!=NULL)
            p[l * NEQ + k] = p[l * NEQ + k] + md[ispec].gir[l*md[ispec].nlev+k];
        }

   //Radiation trapping using the Escape Probability method
   if(par->useEP==1){ 
   double tau, beta, molDens[par->nSpecies];

   molNumDensity(0.0,0.0,subGrid*radius,molDens); 
    for(li=0;li<md[ispec].nline;li++){
      upper=md[ispec].lau[li];
      lower=md[ispec].lal[li];

      //Calculating the optical depth
      tau = ((A[li]*pow(CLIGHT,3))/(8*PI*pow(md[ispec].freq[li],3))) * ((md[ispec].gstat[upper]/md[ispec].gstat[lower])*Pops[lower] - Pops[upper]) * ((molDens[ispec]* radius)/vexp);

      //If the optical depth is small, ignore it
      if (tau > -1.0e-6 && tau < 1.0e-6){
        beta = 1.0;
      }else if (tau>0.0) {
        beta = (2/(3*tau)) - exp(-tau/2)*(tau*(gsl_sf_bessel_Kn(2,tau/2)-gsl_sf_bessel_K1(tau/2))/3 -  gsl_sf_bessel_K1(tau/2)); 
      }else if (tau<0.0){ 
        beta = (1 - exp(-tau)) / tau;
      }

      p[upper * NEQ + lower] = p[upper * NEQ + lower] + A[li]*beta;

    }//end for
  }//end if

  //No photon trapping simulation
  else if(par->useEP==0){
    for(li=0;li<md[ispec].nline;li++){
      upper=md[ispec].lau[li];
      lower=md[ispec].lal[li];
      p[upper * NEQ + lower] = p[upper * NEQ + lower] + A[li];
    }
  
  //Calculate photon trapping self-consistently, using full 3D treatment of radiation field, and photon propagation
  }else if(par->useEP==2){

  if(radius >= gp[gp_sorter[subGrid_pIntensity-1]].radius){
    for(li=0;li<md[ispec].nline;li++){
      upper=md[ispec].lau[li];
      lower=md[ispec].lal[li];
  // Note - this is broken as jbar_grid should refer to the current subgrid, not the total grid, and we should be choosing the closest grid point to the cone center, relative to our present radius     
      p[upper * NEQ + lower] = p[upper * NEQ + lower] + md[ispec].beinstl[li]*jbar_grid[(par->pIntensity-1) * md[ispec].nline + li];
      p[lower * NEQ + upper] = p[lower * NEQ + upper] + md[ispec].beinstl[li]*jbar_grid[(par->pIntensity-1) * md[ispec].nline + li];
    }
  }


  else{
    for(li=0;li<md[ispec].nline;li++){
      upper=md[ispec].lau[li];
      lower=md[ispec].lal[li];
  // Note - this calculation is badly broken as jbar_grid should refer to the current subgrid, not the total grid, and we should be choosing the closest grid point to the cone center, relative to our present radius... (par->pIntensity-1) is used only as a placeholder grid point, to allow the code to run without causing major problems
      p[upper * NEQ + lower] = p[upper * NEQ + lower] + md[ispec].beinstl[li]*jbar_grid[(par->pIntensity-1) * md[ispec].nline + li];
      p[lower * NEQ + upper] = p[lower * NEQ + upper] + md[ispec].beinstl[li]*jbar_grid[(par->pIntensity-1) * md[ispec].nline + li];
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
  if(par->outputfile) popsout(par,gp,md);

}


/*....................................................................*/
/* Sets the Differential Equation (Pdot) to be solved by CVode */
int f(realtype radius, N_Vector P, N_Vector Pdot, void *data){

  int i, j, NEQ, id;
  struct transitionParams *user_data = data;
  NEQ = user_data -> array_size;
  double *p = user_data -> transition_rates;
  double Pops_array[NEQ], vel[DIM], vexp;

  velocity(0.,0.,user_data->subGrid*radius,vel);
  vexp = sqrt(vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2]);

  //Pasing P values to Pops_array for readability
  for(i=0; i < NEQ; ++i)
    Pops_array[i] = Ith(P,i);

  getTransitionRates(user_data->md,user_data->ispec,user_data->gp,user_data->par, NEQ, user_data -> A_array, p, radius, user_data->jbar_grid, Pops_array, user_data->nMaserWarnings, user_data->gp_sorter, user_data->subGrid, user_data->gp_pIntensity, vexp);

  //Initializing Pdot (otherwise it stores previous values between calls to CVode)
  for(i=0; i < NEQ; ++i)
    Ith(Pdot,i) = 0.0;

  //Setting Pdot
  for(i=0; i < NEQ; ++i)
    for(j=0; j < NEQ; ++j){
      if(i != j) 
        Ith(Pdot,i) = Ith(Pdot,i) + (Pops_array[j] * (p[j* NEQ +i])) - (Pops_array[i] * p[i * NEQ +j]);
    }
   
     // Multiply equations by 1/V to get dP/dr     
  for(i=0; i < NEQ; ++i){
     Ith(Pdot,i) = Ith(Pdot,i) / vexp;
    //  Prind the Pdots for debugging
    //  printf("%d %le\n",i,  Ith(Pdot,i));
   }
 
    
  return(0);
}

void reduceTol(N_Vector abstol, realtype reltol, void *cvode_mem, double factor, int NEQ){
   printf("INFO: Reducing RTOL and ATOL by 0.1 and restarting CVODE\n");
   int i, retval;

   for(i=0; i < NEQ; ++i){
      Ith(abstol,i) = Ith(abstol,i) * factor;
  }

   reltol = reltol * factor;
   
  retval = CVodeSVtolerances(cvode_mem, reltol, abstol);
  if (check_retval(&retval, "CVodeSVtolerances", 1)) return;

}


/*....................................................................*/
void
solveStatEq(struct grid *gp, molData *md, const int ispec, configInfo *par\
  , struct blendInfo blends, int *nextMolWithBlend, gridPointData **mp\
  , double **halfFirstDs, int *nMaserWarnings,struct grid *gp3D,int gp_pIntensity, int gp_ncell, double *radii, int subGrid_pIntensity, int *gp_sorter, int subGrid){

  int id;
  realtype reltol, t;
  N_Vector P, abstol;
  SUNMatrix sunMatrix;
  SUNLinearSolver LS;
  void *cvode_mem;
  int i,j,k, cvodeErrs = 0;;
  int retval,cvstatus,index;
  int ncell = par->ncell, pIntensity = par->pIntensity;

  gsl_vector *newpop = gsl_vector_alloc(md[ispec].nlev);

  /* Initializing parameters to be used by CVode */
  int NEQ = md[ispec].nlev; // Number of Equations .
  double p[NEQ][NEQ]; //Collisional rates
  double A[md[ispec].nline]; //Einstein As
  double Pops[NEQ]; //Level Populations
  double popGrid[NRADS][NEQ]; //Populations as a function of radius 
  
  double (*jbar_grid)[md[ispec].nline] = malloc(sizeof(double[pIntensity][md[ispec].nline])); 
   
  for(id=0;id<md[ispec].nline;id++)
    A[id] = md[ispec].aeinst[id];
  
  if (par->useEP==2){
  //Changing the values to those of the full grid to perform the jbar update
  par->ncell = gp_ncell;
  par->pIntensity = gp_pIntensity; 

    for(id=0;id<par->pIntensity;id++){
      updateJBar(id,md,gp3D,ispec,par,blends,nextMolWithBlend[id],mp[id],halfFirstDs[id]);
    }
    for(i=0;i<par->pIntensity;i++){
      for(j=0;j<md[ispec].nline;j++){
          for(k=0;k<gp_pIntensity;k++){
            if(gp3D[i].id==gp[k].id) //we use the .id parameter to to the mapping between the grid points in the full grid and the subgrids
              break;
          }
        jbar_grid[k][j] = mp[i][ispec].jbar[j];
      }
    }
  //Reverting to the original values in case they were changed when useEP==2
  par->ncell = ncell;
  par->pIntensity = pIntensity;
  }

  /*Initializing Pops */
  for(i=0;i<md[ispec].nlev;i++){
    Pops[i] = gp[gp_sorter[0]].mol[ispec].pops[i]; 
    popGrid[0][i] = Pops[i];
  }

  struct transitionParams user_data = {NEQ, A, *p, md,ispec,gp,par, *jbar_grid, nMaserWarnings, gp_sorter, subGrid, subGrid_pIntensity}; 

  P = abstol = NULL;
  sunMatrix = NULL;
  LS = NULL;
  cvode_mem = NULL;

  /* Create serial vector of length NEQ for I.C. and abstol */
  P = N_VNew_Serial(NEQ);
  if (check_retval((void *)P, "N_VNew_Serial", 0)) return;
  abstol = N_VNew_Serial(NEQ); 
  if (check_retval((void *)abstol, "N_VNew_Serial", 0)) return;

  /* Set the scalar relative tolerance */
  reltol = RTOL;

  /* Set the initial population values and the tolerances*/
  for(i=0; i < NEQ; ++i){
    Ith(P,i) = Pops[i];
    Ith(abstol,i) = ATOL;
  }

  /* Call CVodeCreate to create the solver memory and specify the 
  * Backward Differentiation Formula */
  cvode_mem = CVodeCreate(CV_BDF);
  if (check_retval((void *)cvode_mem, "CVodeCreate", 0)) return;

  /* Call CVodeInit to initialize the integrator memory and specify the
   * user's right hand side function in y'=f(t,y), the inital time T0, and
   * the initial dependent variable vector y. */
  retval = CVodeInit(cvode_mem, f,  radii[0], P);
  if (check_retval(&retval, "CVodeInit", 1)) return;

  retval = CVodeSetUserData(cvode_mem, &user_data);
  if (check_retval(&retval, "CVodeSetUserData", 0)) return;

  /* Call CVodeSVtolerances to specify the scalar relative tolerance
   * and vector absolute tolerances */
  retval = CVodeSVtolerances(cvode_mem, reltol, abstol);
  if (check_retval(&retval, "CVodeSVtolerances", 1)) return;

  retval = CVodeSetMaxNumSteps(cvode_mem, 5000);
    if (check_retval(&retval, "CVodeSetMaxNumSteps", 1)) return;

  /* Create dense SUNMatrix for use in linear solves */
  sunMatrix = SUNDenseMatrix(NEQ, NEQ);
  if(check_retval((void *)sunMatrix, "SUNDenseMatrix", 0)) return;

  /* Create dense SUNLinearSolver object for use by CVode */
  LS = SUNLinSol_Dense(P, sunMatrix);
  if(check_retval((void *)LS, "SUNLinSol_Dense", 0)) return;

  /* Call CVodeSetLinearSolver to attach the matrix and linear solver to CVode */
  retval = CVodeSetLinearSolver(cvode_mem, LS, sunMatrix);
  if(check_retval(&retval, "CVodeSetLinearSolver", 1)) return;

  /* Set the user-supplied Jacobian routine Jac */
  //retval = CVodeSetJacFn(cvode_mem, Jac);
  //if(check_retval(&retval, "CVodeSetJacFn", 1)) return;

  /* Use a difference quotient Jacobian */
  retval = CVodeSetJacFn(cvode_mem, NULL);
  if(check_retval(&retval, "CVodeSetJacFn", 1)) return;
  
  printf("Starting CVODE time-dependent solver for sub-grid %d, species %d...\n",subGrid,ispec);
  fflush(stdout);

  // Call CVODE for each radius
  // Do/while loop to catch CVODE error test failure status
  do{
     for(i=1; i<NRADS; i++){
         cvstatus = CVode(cvode_mem, radii[i], P, &t, CV_NORMAL);
        
       if(cvstatus == CV_SUCCESS){
         for(j=0;j<md[ispec].nlev;j++){ 
           popGrid[i][j] = Ith(P,j); 
         }
       }
      
      if(cvstatus==-3 && reltol > MINTOL){ // Reduce the tolerances and try again
         retval = CVodeInit(cvode_mem, f,  radii[0], P);
         if (check_retval(&retval, "CVodeInit", 1)) printf("Failed to reinitialize CVODE\n");
         CVodeSetUserData(cvode_mem, &user_data);
         CVodeSetMaxNumSteps(cvode_mem, 5000);
         sunMatrix = SUNDenseMatrix(NEQ, NEQ);
         LS = SUNLinSol_Dense(P, sunMatrix);
         CVodeSetLinearSolver(cvode_mem, LS, sunMatrix);
         CVodeSetJacFn(cvode_mem, NULL);
         reduceTol(abstol, reltol, cvode_mem, 0.1, NEQ);
         break;
      }else if(cvstatus!=0){
      cvodeErrs++;
      printf("CVODE error %d (continuing to next timestep - check populations!!)\n",cvstatus);
      }
    //Next CVODE step
      if(cvodeErrs>=15){ 
         bail_out("CVODE solver failure - check physical model and tolerances.");
         exit(1);
      }
    }
  }while(cvstatus==-3) ;

  fflush(stdout);
  /*Interpolate the computed radial populations onto the Delaunay grid for raytracing*/ 
  index=0;
  for(i=0;i<subGrid_pIntensity;i++){
     for(j=1;j<NRADS;j++){ 
        if(radii[j]>gp[i].radius){
           index = j;
           break;
        }
     }   
    fflush(stdout);     
     for(k=0;k<md[ispec].nlev;k++){
        gp[i].mol[ispec].pops[k] = linear_interp(radii[index-1],radii[index],popGrid[index-1][k],popGrid[index][k],gp[i].radius); 
     }
  }

  /* Free P and abstol vectors */
  N_VDestroy(P);
  N_VDestroy(abstol);

  /* Free integrator memory */
  CVodeFree(&cvode_mem);

  /* Free the linear solver memory */
  SUNLinSolFree(LS);

  /* Free the matrix memory */
  SUNMatDestroy(sunMatrix);
  
  gsl_vector_free(newpop);
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
   if(par->outputfile && par->useEP != 2) popsout(par,gp,md);

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

    for(id=0;id<par->pIntensity;id++){
      totalNMaserWarnings = nMaserWarnings[id];
    }

    if(!silent && totalNMaserWarnings>0){
      snprintf(message, STR_LEN_0, "Maser warning: optical depth dropped below -%4.1f %d times this iteration.", MAX_NEG_OPT_DEPTH, totalNMaserWarnings);
      warning(message);
    }

    if(!silent) warning("");
    for (i=0;i<par->pIntensity;i++){
      freeGridPointData(par->nSpecies, mp[i]);
      free(halfFirstDs[i]);
    }
    if(par->outputfile != NULL) popsout(par,gp,md);

  if(par->useEP == 2) freeMolsWithBlends(blends.mols, blends.numMolsWithBlends);
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