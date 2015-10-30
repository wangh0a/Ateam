#include "sinc.h"
#ifdef _OPENMP
#include <omp.h>
#endif

/*
 * Possible references:
 *   Arbitrary source and receiver positioning in finite-difference 
 *   schemes using Kaiser windowed sinc functions
 *   2002, Hicks
 * Notice:
 *    9-point interpolation along each axis
 *    source and receiver position must be 4 grids away from the boundary
 */


/*------------------------------------------------------------*/
scoef2d sinc2d_make(int nc,
    pt2d* aa,
    int nz, int nx,
    float dz, float dx,
    float oz, float ox)
/*< init the sinc2d interpolation for injection/extraction >*/
/* 8-point neighbors on each axis */
{
  int i, ic;
  scoef2d swout;

  float inp[9]; 
  float xo[9];

  int ix, iz;
  float distx, distz;

  swout = (scoef2d) sf_alloc(nc,sizeof(*swout));

  for (i=0; i<9; i++) inp[i] = 0.0f;
  inp[4] = 1.0f;
  /* allocate and set loop */
  for (ic=0; ic<nc; ++ic){
    ix = (int)((aa[ic].x - ox)/dx+0.499f); 
    iz = (int)((aa[ic].z - oz)/dz+0.499f);
    swout[ic].fx = 0;
    swout[ic].nx = 9;
    swout[ic].fz = 0;
    swout[ic].nz = 9;
    swout[ic].n = nc;

    swout[ic].ix = ix;
    swout[ic].iz = iz;

    distx = ix*dx + ox - aa[ic].x; 
    distz = iz*dz + oz - aa[ic].z;

    for (i=0; i<9; ++i)  xo[i] = - 4.0f + distx/dx + i*1.0f;
    ints8r (9, 1.0f, -4.0, inp, 0.0f, 0.0f, 9, xo, swout[ic].sincx);

    if(swout[ic].sincx[4]==1.0f){
      swout[ic].nx = 1;
      swout[ic].fx = 4;
    }

    for (i=0; i<9; ++i)  xo[i] = - 4.0 + distz/dz + i*1.0f;
    ints8r (9, 1.0f, -4.0, inp, 0.0f, 0.0f, 9, xo, swout[ic].sincz);

    if(swout[ic].sincz[4]==1.0f){
      swout[ic].nz = 1;
      swout[ic].fz = 4;
    }

  }  
  return swout;
}


/*------------------------------------------------------------*/
void sinc2d_inject(float**uu,
    float *dd,
    scoef2d ca)
/*< inject into wavefield >*/
{

  int   ia, ix, iz, sx, sz, ixx, izz;
  float w, wx, wz;
  float value;

  int na = ca[0].n;
#ifdef _OPENMP
#pragma omp parallel for      \
  schedule(dynamic)       \
  private(ia,ix,iz,sx,sz,w,wx,wz,ixx,izz,value)   \
  shared(na,ca,dd,uu)
#endif
  for(ia=0;ia<na;ia++) {  
    w = dd[ia];
    ix = ca[ia].ix;
    iz = ca[ia].iz;
    for (ixx=ca[ia].fx; ixx<ca[ia].fx+ca[ia].nx; ixx++){
      sx = -4 +ixx;
      wx = ca[ia].sincx[ixx];
      for(izz=ca[ia].fz; izz<ca[ia].fz+ca[ia].nz; izz++){
        sz = -4 +izz;
        wz = ca[ia].sincz[izz];
        value = w*wx*wz;
#ifdef _OPENMP
#pragma omp atomic
#endif
        uu[ix+sx][iz+sz] += value; /* scatter */
      }
    }  
  }
}

/*------------------------------------------------------------*/
void sinc2d_inject1(float**uu,
    float dd,
    scoef2d ca)
/*< inject into wavefield >*/
{

  int   ia, ix, iz, sx, sz, ixx, izz;
  float w, wx, wz;
  float value;
  int na = ca[0].n;

#ifdef _OPENMP
#pragma omp parallel for      \
  schedule(dynamic)       \
  private(ia,ix,iz,sx,sz,w,wx,wz,ixx,izz,value)   \
  shared(na,ca,dd,uu)
#endif
  for(ia=0;ia<na;ia++) {
    w = dd;
    ix = ca[ia].ix;
    iz = ca[ia].iz;
    for (ixx=ca[ia].fx; ixx<ca[ia].fx+ca[ia].nx;ixx++){
      sx = -4 +ixx;
      wx = ca[ia].sincx[ixx];
      for(izz=ca[ia].fz; izz<ca[ia].fz+ca[ia].nz; izz++){
        sz = -4 +izz;
        wz = ca[ia].sincz[izz];
        value = w*wx*wz;
#ifdef _OPENMP
#pragma omp atomic
#endif
        uu[ix+sx][iz+sz] += value; /* scatter */
      }
    }  
  }
}

/*------------------------------------------------------------*/
void sinc2d_extract(float**uu,
    float *dd,
    scoef2d ca)
/*< inject into wavefield >*/
{
  int   ia, ix, iz, sx, sz, ixx, izz;
  float wx, wz;
  float gather;
  int na = ca[0].n;

#ifdef _OPENMP
#pragma omp parallel for      \
  schedule(dynamic)       \
  private(ia,ix,iz,sx,sz,wx,wz,ixx,izz,gather)    \
  shared(na,ca,dd,uu)
#endif
  for(ia=0;ia<na;ia++) {
    ix = ca[ia].ix;
    iz = ca[ia].iz;
    gather = 0.f;
    for (ixx=ca[ia].fx; ixx<ca[ia].fx+ca[ia].nx; ixx++){
      sx = -4 +ixx;
      wx = ca[ia].sincx[ixx];
      for(izz=ca[ia].fz; izz<ca[ia].fz+ca[ia].nz; izz++){
        sz = -4 +izz;
        wz = ca[ia].sincz[izz];
        gather += uu[ix+sx][iz+sz]*wx*wz; /* gather */
      }
    }  
    dd[ia] = gather;
  }
}


/* == 3D == */
/*------------------------------------------------------------*/
scoef3d sinc3d_make(int nc,
    pt3d* aa,
    int nz, int nx, int ny,
    float dz, float dx, float dy,
    float oz, float ox, float oy)
/*< init the sinc3d interpolation for injection/extraction >*/
/* 8-point neighbors on each axis */
{
  int i, ic;
  scoef3d swout;

  float inp[9]; 
  float xo[9];

  int iy, ix, iz;
  float disty, distx, distz;

  swout = (scoef3d) sf_alloc(nc,sizeof(*swout));

  for (i=0; i<9; i++) inp[i] = 0.0f;
  inp[4] = 1.0f;
  /* allocate and set loop */
  for (ic=0; ic<nc; ++ic){
    iy = (int)((aa[ic].y - oy)/dy+0.499f); 
    ix = (int)((aa[ic].x - ox)/dx+0.499f); 
    iz = (int)((aa[ic].z - oz)/dz+0.499f);
    swout[ic].fy = 0;
    swout[ic].ny = 9;
    swout[ic].fx = 0;
    swout[ic].nx = 9;
    swout[ic].fz = 0;
    swout[ic].nz = 9;
    swout[ic].n = nc;

    swout[ic].iy = iy;
    swout[ic].ix = ix;
    swout[ic].iz = iz;

    disty = iy*dy + oy - aa[ic].y; 
    distx = ix*dx + ox - aa[ic].x; 
    distz = iz*dz + oz - aa[ic].z;

    /* y direction */
    for (i=0; i<9; ++i)  xo[i] = - 4.0f + disty/dy + i*1.0f;
    ints8r (9, 1.0f, -4.0, inp, 0.0f, 0.0f, 9, xo, swout[ic].sincy);
    if(swout[ic].sincy[4]==1.0f){
      swout[ic].ny = 1;
      swout[ic].fy = 4;
    }
    /* x direction */
    for (i=0; i<9; ++i)  xo[i] = - 4.0f + distx/dx + i*1.0f;
    ints8r (9, 1.0f, -4.0, inp, 0.0f, 0.0f, 9, xo, swout[ic].sincx);

    if(swout[ic].sincx[4]==1.0f){
      swout[ic].nx = 1;
      swout[ic].fx = 4;
    }

    /* z direction */
    for (i=0; i<9; ++i)  xo[i] = - 4.0 + distz/dz + i*1.0f;
    ints8r (9, 1.0f, -4.0, inp, 0.0f, 0.0f, 9, xo, swout[ic].sincz);

    if(swout[ic].sincz[4]==1.0f){
      swout[ic].nz = 1;
      swout[ic].fz = 4;
    }
  }  
  return swout;
}


/*------------------------------------------------------------*/
void sinc3d_inject(float ***uu,
    float *dd,
    scoef3d ca)
/*< inject into wavefield >*/
{

  int   ia, iy, ix, iz, sy, sx, sz, iyy, ixx, izz;
  float w, wy, wx, wz;
  float value;

  int na = ca[0].n;
#ifdef _OPENMP
#pragma omp parallel for      \
  schedule(dynamic)       \
  private(ia,iy,ix,iz,sy,sx,sz,w,wy,wx,wz,iyy,ixx,izz,value)   \
  shared(na,ca,dd,uu)
#endif
  for(ia=0;ia<na;ia++) {  
    w = dd[ia];
    iy = ca[ia].iy;
    ix = ca[ia].ix;
    iz = ca[ia].iz;
    for (iyy=ca[ia].fy; iyy<ca[ia].fy+ca[ia].ny; iyy++){
      sy = -4 +iyy;
      wy = ca[ia].sincy[iyy];
      for (ixx=ca[ia].fx; ixx<ca[ia].fx+ca[ia].nx; ixx++){
        sx = -4 +ixx;
        wx = ca[ia].sincx[ixx];
        for(izz=ca[ia].fz; izz<ca[ia].fz+ca[ia].nz; izz++){
          sz = -4 +izz;
          wz = ca[ia].sincz[izz];
          value = w*wy*wx*wz;
#ifdef _OPENMP
#pragma omp atomic
#endif
          uu[iy+sy][ix+sx][iz+sz] += value; /* scatter */
        }
      }
    }
  }
}

/*------------------------------------------------------------*/
void sinc3d_inject1(float ***uu,
    float dd,
    scoef3d ca)
/*< inject into wavefield >*/
{

  int   ia, iy, ix, iz, sy, sx, sz, iyy, ixx, izz;
  float w, wy, wx, wz;
  float value;
  int na = ca[0].n;

#ifdef _OPENMP
#pragma omp parallel for      \
  schedule(dynamic)       \
  private(ia,iy,ix,iz,sy,sx,sz,w,wy,wx,wz,iyy,ixx,izz,value)   \
  shared(na,ca,dd,uu)
#endif
  for(ia=0;ia<na;ia++) {
    w = dd;
    iy = ca[ia].iy;
    ix = ca[ia].ix;
    iz = ca[ia].iz;
    for (iyy=ca[ia].fy;iyy<ca[ia].fy+ca[ia].ny;iyy++){
      sy = -4 +iyy;
      wy = ca[ia].sincy[iyy];
      for (ixx=ca[ia].fx; ixx<ca[ia].fx+ca[ia].nx;ixx++){
        sx = -4 +ixx;
        wx = ca[ia].sincx[ixx];
        for(izz=ca[ia].fz; izz<ca[ia].fz+ca[ia].nz; izz++){
          sz = -4 +izz;
          wz = ca[ia].sincz[izz];
          value = w*wy*wx*wz;
#ifdef _OPENMP
#pragma omp atomic
#endif
          uu[iy+sy][ix+sx][iz+sz] += value; /* scatter */
        }
      }
    }
  }
}

/*------------------------------------------------------------*/
void sinc3d_extract(float ***uu,
    float *dd,
    scoef3d ca)
/*< inject into wavefield >*/
{
  int   ia, iy, ix, iz, sy, sx, sz, iyy, ixx, izz;
  float wy, wx, wz;
  float gather;
  int na = ca[0].n;

#ifdef _OPENMP
#pragma omp parallel for      \
  schedule(dynamic)       \
  private(ia,iy,ix,iz,sy,sx,sz,wy,wx,wz,iyy,ixx,izz,gather)    \
  shared(na,ca,dd,uu)
#endif
  for(ia=0;ia<na;ia++) {
    iy = ca[ia].iy;
    ix = ca[ia].ix;
    iz = ca[ia].iz;
    gather = 0.f;
    for (iyy=ca[ia].fy; iyy<ca[ia].fy+ca[ia].ny; iyy++){
      sy = -4 +iyy;
      wy = ca[ia].sincy[iyy];
      for (ixx=ca[ia].fx; ixx<ca[ia].fx+ca[ia].nx; ixx++){
        sx = -4 +ixx;
        wx = ca[ia].sincx[ixx];
        for(izz=ca[ia].fz; izz<ca[ia].fz+ca[ia].nz; izz++){
          sz = -4 +izz;
          wz = ca[ia].sincz[izz];
          gather += uu[iy+sy][ix+sx][iz+sz]*wy*wx*wz; /* gather */
        }
      }
    }
    dd[ia] = gather;
  }
}



