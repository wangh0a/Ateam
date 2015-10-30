#ifndef _SINC_H
#define _SINC_H

#include "rsf.h"
#include "rsf_su.h"

typedef struct scoef2 *scoef2d;
struct scoef2{
  int n;
  int ix,iz;
  int fx,fz,nx,nz;
  float sincx[9], sincz[9];
};

typedef struct scoef3 *scoef3d;
struct scoef3{
  int n;
  int iy,ix,iz;
  int fy,fx,fz,ny,nx,nz;
  float sincy[9], sincx[9], sincz[9];
};

scoef2d sinc2d_make(int nc, pt2d *aa, int nz, int nx, float dz, float dx, float oz, float ox);
void sinc2d_inject(float **uu, float *dd, scoef2d ca);
void sinc2d_inject1(float **uu, float dd, scoef2d ca);
void sinc2d_extract(float **uu, float *dd, scoef2d ca);

scoef3d sinc3d_make(int nc, pt3d *aa, int nz, int nx, int ny, float dz, float dx, float dy, float oz, float ox, float oy);
void sinc3d_inject(float ***uu, float *dd, scoef3d ca);
void sinc3d_inject1(float ***uu, float dd, scoef3d ca);
void sinc3d_extract(float ***uu, float *dd, scoef3d ca);
#endif
