#include "util.h"

/* --------------------------------------------------- */
void
expand2d(float *uo, float *ui, int nzpad, int nxpad, int nz, int nx)
{
  int nb = 0.5 * (nzpad - nz);
  for (int ix=0; ix<nx; ix++) {
    for (int iz=0; iz<nz; iz++) {
      uo[(ix+nb)*nzpad+iz+nb] = ui[ix*nz+iz];
    }
  }

  for (int ib=0; ib<nb; ib++) {
    for (int iz=0; iz<nzpad; iz++) {
      uo[ib*nzpad+iz] = uo[nb*nzpad+iz];
      uo[(nx+nb+ib)*nzpad+iz] = uo[(nx+nb-1)*nzpad+iz];
    }
  }

  for (int ix=0; ix<nxpad; ix++) {
    for (int ib=0; ib<nb; ib++) {
      uo[ix*nzpad+ib] = uo[ix*nzpad+nb];
      uo[ix*nzpad+nz+nb+ib] = uo[ix*nzpad+nz+nb-1];
    }
  }
  return;
}

/* --------------------------------------------------- */
void
wwin2d(float **uo, float **ui, int nzo, int nxo, int nb)
{
#pragma omp parallel for schedule(dynamic, 1) 
  for (int ix=0; ix<nxo; ix++) {
    for (int iz=0; iz<nzo; iz++) {
      uo[ix][iz] = ui[ix+nb][iz+nb];
    }
  }
  return;
}

/* --------------------------------------------------- */
void
wfld2d_inject(float **uo, float **ui, int nzo, int nxo, int nb)
{
#pragma omp parallel for schedule(dynamic, 1) 
  for (int ix=0; ix<nxo; ix++) {
    for (int iz=0; iz<nzo; iz++) {
      uo[ix+nb][iz+nb] += ui[ix][iz];
    }
  }
  return;
}


/* --------------------------------------------------- */
float *damp_make(int nb)
{
  float *damp = NULL;
  // if (nb > 0) damp = sf_floatalloc(nb);
  if (nb > 0) damp = malloc(nb*sizeof(float));
  /* if (nb > 0) damp = new float[nb]; */
  for (int i=0; i<nb; i++) {
    float fb = 1.f * i / nb;
    damp[nb-1-i] = expf(- fb * fb / 8.f);
    // sf_warning("i=%d damp=%f", i, damp[nb-1-i]);
  }
  return damp;
}

/* --------------------------------------------------- */
void
damp2d_apply(float **uu, float *damp, int nz, int nx, int nb)
{
  if (damp != NULL) {
#pragma omp parallel for schedule(dynamic, 1) 
    for (int ib=0; ib<nb; ib++) {
      /* top/bottom boundary */
      for (int iz=0; iz<nz; iz++) {
        uu[ib][iz] *= damp[ib];
        uu[nx-1-ib][iz] *= damp[ib];
      }
    }
#pragma omp parallel for schedule(dynamic, 1) 
    for (int ix=0; ix<nx; ix++) {
      for (int ib=0; ib<nb; ib++) {
        /* left/right boundary */
        uu[ix][ib] *= damp[ib];
        uu[ix][nz-1-ib] *= damp[ib];
      }
    }
  }
  return;
}

