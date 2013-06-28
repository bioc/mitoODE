#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rdefines.h>
#include <R_ext/Error.h>

void rksolve(int nt, double *vt, int nc, double *vy0,
             void (*gf)(double t, double *vy, double avy, double *tvy, double *vp, double *gvy),
	     double *vp, double *my);

#define MAX_ODEGF 5
void gf0(double t, double *vy, double avy, double *tvy, double *vp, double *gvy);
void gf1(double t, double *vy, double avy, double *tvy, double *vp, double *gvy);
void gf2(double t, double *vy, double avy, double *tvy, double *vp, double *gvy);
void gf3(double t, double *vy, double avy, double *tvy, double *vp, double *gvy);
void gf4(double t, double *vy, double avy, double *tvy, double *vp, double *gvy);
void *ODEgf[MAX_ODEGF] = { &gf0, &gf1, &gf2, &gf3, &gf4};

/**
 *  ODE Runge-Kutta 4th order solver
 *
 *  input:
 *  vt: vector of time points (of size nt)
 *  vy0: vector of initial values (of size nc)
 *  vp: vector of ODE gradient parameters (of size np)
 *  igf: integer indexing which ODE gradient function among ODEgf should be used
 *
 *  output:
 *  return: solution matrix (of size nt*nc)
 *
 *  reference:
 *  http://mathworld.wolfram.com/Runge-KuttaMethod.html
 */
SEXP rksolve_wrap(SEXP _vt, SEXP _vy0, SEXP _vp, SEXP _igf) {
  int nt,nc;
  double *vt;
  double *vy0;
  double *vp;
  int igf;
  double *my;
  int nprotect=0;
  SEXP result, resultdim;

  nt = LENGTH(_vt);
  vt = REAL(_vt);

  nc = LENGTH(_vy0);
  vy0 = REAL(_vy0);
  vp = REAL(_vp);
  igf = INTEGER(_igf)[0];
  if (igf<0 || igf>=MAX_ODEGF) error( "'igf' must index a known ODE gradient functions (see rksolve.c)");

  PROTECT( result=allocVector(REALSXP, nt*nc) );
  nprotect++;
  PROTECT( resultdim=allocVector(INTSXP, 2) );
  nprotect++;

  my=REAL(result);
  INTEGER(resultdim)[0]=nt;
  INTEGER(resultdim)[1]=nc;
  SET_DIM(result, resultdim);
  
  rksolve(nt, vt, nc, vy0, ODEgf[igf], vp, my);

  UNPROTECT(nprotect);
  return result;
}

/**
 *  ODE Runge-Kutta 4th order solver
 *
 *  input:
 *  nt: number of time points
 *  vt: vector of time points (of size nt)
 *  nc: number of components of the solution
 *  vy0: vector of initial values (of size nc)
 *  gf: ODE gradient function
 *  vp: vector of ODE gradient parameters (of size np)
 *
 *  output:
 *  my: matrix solution (must be allocated before call, of size nt*nc)
 *
 *  reference:
 *  http://mathworld.wolfram.com/Runge-KuttaMethod.html
 */
void rksolve(int nt, double *vt, int nc, double *vy0, 
	     void (*gf)(double t, double *vy, double avy, double *tvy, double *vp, double *gvy),
	     double *vp, double *my) {
  int i, j;
  double t, dt;
  double *mk;
  double *vy;

  mk=(double *)malloc(sizeof(double)*nc*4);
  vy=(double *)malloc(sizeof(double)*nc);

  // copy initial values
  for (j=0;j<nc;j++) vy[j]=vy0[j];
  for (j=0;j<nc;j++) my[0+j*nt]=vy0[j];

  // runge-kutta 4th order
  for (i=1;i<nt;i++) {
    t=vt[i];
    dt=vt[i]-vt[i-1];

    gf(t     ,vy,0.0   ,NULL     ,vp,&mk[0]);
    gf(t+dt/2,vy,dt*0.5,&mk[0]   ,vp,&mk[nc]);
    gf(t+dt/2,vy,dt*0.5,&mk[nc]  ,vp,&mk[2*nc]);
    gf(t+dt  ,vy,dt    ,&mk[2*nc],vp,&mk[3*nc]);
    
    for (j=0;j<nc;j++) vy[j]+=dt*(mk[j]+2*mk[j+nc]+2*mk[j+2*nc]+mk[j+3*nc])/6;
    for (j=0;j<nc;j++) my[i+j*nt]=vy[j];  
  }

  free(vy);
  free(mk);
}

/**
 *  ODE gradient function, gf(vy+tvy*avy, t, vp)
 * 
 *  input:
 *  vy: vector y
 *  t: time point
 *  vp: vector of parameters
 *  tvy, avy: if tvy non-null, derivate will be evaluated at vy + avy*tvy instead of vy 
 *
 *  output:
 *  gvy: gradient (must be allocated before)
 */
void gf0(double t, double *vy, double avy, double *tvy, double *vp, double *gvy) {
  double ni=vy[0];
  double nm=vy[1];
  double ns=vy[2];
  double na=vy[3];
  double u =vy[4];

  double k1 = vp[0];
  double k2 = vp[1];
  double k3 = vp[2];
  double k4 = vp[3];
  double k5 = vp[4];
  double k6 = vp[5];
  double k7 = vp[6];

  // Change the valuation point
  if (tvy!=NULL) {
    ni+=avy*tvy[0];
    nm+=avy*tvy[1];
    ns+=avy*tvy[2];
    na+=avy*tvy[3];
    u +=avy*tvy[4];
  }
  
  gvy[0] = -k1*u*ni -k2*ni + 2*k3*nm;
  gvy[1] = k1*u*ni - (k3+k4+k5)*nm;
  gvy[2] = k5*nm - k6*ns;
  gvy[3] = k2*ni + k4*nm + k6*ns;
  gvy[4] = -k7*u*ni;
}

void gf1(double t, double *vy, double avy, double *tvy, double *vp, double *gvy) {
  double ni = vy[0];
  double nm = vy[1];
  double ns = vy[2];
  double na = vy[3];

  double k1 = vp[0] / (1 + exp(0.3*(t-vp[1])));
  double k2 = vp[2] / (1 + exp(-0.3*(t-vp[3])));
  double k3 = vp[4] / (1 + exp(0.3*(t-vp[5])));
  double k4 = vp[6] / (1 + exp(-0.3*(t-vp[7])));
  double k5 = vp[8] / (1 + exp(-0.3*(t-vp[9])));
  double k6 = vp[10];

  // Change the valuation point
  if (tvy!=NULL) {
    ni += avy*tvy[0];
    nm += avy*tvy[1];
    ns += avy*tvy[2];
    na += avy*tvy[3];
  }
  
  gvy[0] = -k1*ni -k2*ni + 2*k3*nm;
  gvy[1] = k1*ni - (k3+k4+k5)*nm;
  gvy[2] = k5*nm - k6*ns;
  gvy[3] = k2*ni + k4*nm + k6*ns;
}

void gf2(double t, double *vy, double avy, double *tvy, double *vp, double *gvy) {
  double ni = vy[0];
  double nm = vy[1];
  double ns = vy[2];
  double na = vy[3];

  double k1 = vp[0] / (1 + exp(0.3*(t-vp[1])));
  double k2 = vp[2] / (1 + exp(-0.3*(t-vp[3])));
  double k3 = vp[4] / (1 + exp(0.3*(t-vp[5])));
  double k4 = vp[6] / (1 + exp(-0.3*(t-vp[7])));
  double k5 = vp[8] / (1 + exp(-0.3*(t-vp[9])));
  double k6 = vp[10] / (1 + exp(-0.3*(t-vp[11])));

  // Change the valuation point
  if (tvy!=NULL) {
    ni += avy*tvy[0];
    nm += avy*tvy[1];
    ns += avy*tvy[2];
    na += avy*tvy[3];
  }
  
  gvy[0] = -k1*ni -k2*ni + 2*k3*nm;
  gvy[1] = k1*ni - (k3+k4+k5)*nm;
  gvy[2] = k5*nm - k6*ns;
  gvy[3] = k2*ni + k4*nm + k6*ns;
}

void gf3(double t, double *vy, double avy, double *tvy, double *vp, double *gvy) {
  double ni = vy[0];
  double nm = vy[1];
  double np = vy[2];
  double na = vy[3];
  const double gamma=1;

  // 0him, hmi, hia, hmp, hma, 5hpa ;  6tim, tmi, ta, tmp, mu, i0, 12aim, ami
  double kim = vp[12] + vp[0] / (1 + exp(-gamma*(t-vp[6])));
  double kmi = vp[13] + vp[1] / (1 + exp(-gamma*(t-vp[7])));
  double kia = vp[2] / (1 + exp(-gamma*(t-vp[8])));
  double kmp = vp[3] / (1 + exp(-gamma*(t-vp[9])));
  double kma = vp[4] / (1 + exp(-gamma*(t-vp[8])));
  double kpa = vp[5] / (1 + exp(-gamma*(t-vp[8])));

  // Change the valuation point
  if (tvy!=NULL) {
    ni += avy*tvy[0];
    nm += avy*tvy[1];
    np += avy*tvy[2];
    na += avy*tvy[3];
  }
  
  gvy[0] = -(kim+kia)*ni + 2*kmi*nm;
  gvy[1] = -(kmi+kmp+kma)*nm + kim*ni;
  gvy[2] = -kpa*np + kmp*nm;
  gvy[3] = kpa*np + kma*nm + kia*ni;
}

void gf4(double t, double *vy, double avy, double *tvy, double *vp, double *gvy) {
  double ni = vy[0];
  double nm = vy[1];
  double np = vy[2];
  double na = vy[3];
  const double gamma=1;

  // 0him, hmi, hmp, ha ; 4tim, tmi, tmp, ta ; 8mu, i0 ; 10aim, 11ami
  double kim = fabs(vp[10] + vp[0] / (1 + exp(-gamma*(t-vp[4]))));
  double kmi = fabs(vp[11] + vp[1] / (1 + exp(-gamma*(t-vp[5]))));
  double kmp = vp[2] / (1 + exp(-gamma*(t-vp[6])));
  double ka  = vp[3] / (1 + exp(-gamma*(t-vp[7])));

  // Change the valuation point
  if (tvy!=NULL) {
    ni += avy*tvy[0];
    nm += avy*tvy[1];
    np += avy*tvy[2];
    na += avy*tvy[3];
  }
  
  gvy[0] = -(kim+ka)*ni + 2*kmi*nm;
  gvy[1] = -(kmi+kmp+ka)*nm + kim*ni;
  gvy[2] = -ka*np + kmp*nm;
  gvy[3] = ka*np + ka*nm + ka*ni;
}
