// DONE: linear approximation (equation 14), with absorbing boundary
//  and  adjoint mode
#include <rsf.h>
#include <omp.h>
#include <complex.h>
#include <fftw3.h>
#include "util.h"
#include "sinc.h"

#define Q1 .5f
#define Q2 .25f

static void
fft_stepforward(
    float **u0, float **u1,
    float *rwave, float *rwavem,
    fftwf_complex *cwave, fftwf_complex *cwavem,
    float **vp, float **vn , float **eta,
    float **vh, float **eps, float **lin_eta,
    float *kz, float *kx,
    fftwf_plan forward_plan, fftwf_plan inverse_plan,
    int nz, int nx, int nzpad, int nxpad,
    int nkz, int nkx,
    float wt, bool adj);

int main (int argc, char *argv[])
{
  bool verb, snap;
  bool abc, adj;
  int nz, nx, nt, ns, nr;
  float dz, dx, dt, oz, ox;
  int nz0, nx0, nb;
  float oz0, ox0;
  int nkz, nkx;
  int nzpad, nxpad;
  
  float **u1, **u0;
  float *ws, *wr;
  
  sf_file file_src = NULL, file_rec = NULL;
  sf_file file_inp = NULL, file_out = NULL;
  sf_file file_mdl = NULL;
  sf_axis az = NULL, ax = NULL, at = NULL, as = NULL, ar = NULL;
  pt2d *src2d = NULL;
  pt2d *rec2d = NULL;
  scoef2d cssinc = NULL;
  scoef2d crsinc = NULL;
  float *wi = NULL, *wo = NULL;
  sf_axis ai = NULL, ao = NULL;
  scoef2d cisinc = NULL, cosinc = NULL;
  bool spt = false, rpt = false;
  bool ipt = false, opt = false;
  
  sf_init(argc, argv);
  
  if (!sf_getbool("verb", &verb)) verb = false;
  if (!sf_getbool("snap", &snap)) snap = false;
  if (!sf_getbool("adj", &adj)) adj = false;
  if (!sf_getint("nb", &nb)) nb = 4;
  if (sf_getstring("sou") != NULL) { 
    spt = true;
    if (adj) opt = true;
    else     ipt = true;
  }
  if (sf_getstring("rec") != NULL) {
    rpt = true;
    if (adj) ipt = true;
    else     opt = true;
  }
  
  file_inp = sf_input("in");
  file_mdl = sf_input("model");
  if (spt) file_src = sf_input("sou");
  if (rpt) file_rec = sf_input("rec");
  file_out = sf_output("out");

  if (ipt) at = sf_iaxa(file_inp, 2);
  else     at = sf_iaxa(file_inp, 3);
  if (spt) as = sf_iaxa(file_src, 2);
  if (rpt) ar = sf_iaxa(file_rec, 2);
  az = sf_iaxa(file_mdl, 1);
  ax = sf_iaxa(file_mdl, 2);
  nt = sf_n(at);  dt = sf_d(at);  //ot = sf_o(at);
  nz0 = sf_n(az);  dz = sf_d(az);  oz0 = sf_o(az);
  nx0 = sf_n(ax);  dx = sf_d(ax);  ox0 = sf_o(ax);

  if (spt) ns = sf_n(as);
  if (rpt) nr = sf_n(ar);
  nz = nz0 + 2 * nb;
  nx = nx0 + 2 * nb;
  oz = oz0 - nb * dz;
  ox = ox0 - nb * dx;
  abc = nb ? true : false;
  // sf_error("ox=%f ox0=%f oz=%f oz0=%f",ox,ox0,oz,oz0);
  
  nzpad = kiss_fft_next_fast_size( ((nz+1)>>1)<<1 );
  nkx = nxpad = kiss_fft_next_fast_size(nx);
  nkz = nzpad / 2 + 1;
  /* float okx = - 0.5f / dx; */
  float okx = 0.f;
  float okz = 0.f;
  float dkx = 1.f / (nxpad * dx);
  float dkz = 1.f / (nzpad * dz);

  float **vp, **vn, **eta;
  vp  = sf_floatalloc2(nz, nx);
  vn  = sf_floatalloc2(nz, nx);
  eta = sf_floatalloc2(nz, nx);
  float **tmparray = sf_floatalloc2(nz0, nx0);
  sf_floatread(tmparray[0], nz0*nx0, file_mdl); expand2d(vp[0], tmparray[0], nz, nx, nz0, nx0);
  sf_floatread(tmparray[0], nz0*nx0, file_mdl); expand2d(vn[0], tmparray[0], nz, nx, nz0, nx0);
  sf_floatread(tmparray[0], nz0*nx0, file_mdl); expand2d(eta[0], tmparray[0], nz, nx, nz0, nx0);

  float **eps, **lin_eta, **vh;  
  eps = NULL, lin_eta = NULL, vh = NULL;
 
  eps     = sf_floatalloc2(nz, nx);
  lin_eta = sf_floatalloc2(nz, nx);
  vh      = sf_floatalloc2(nz, nx);

  float ONE_P_2_DELTA = 0.0f;

  for (int ix=0; ix<nx; ix++) {
    for (int iz=0; iz<nz; iz++){
      vp[ix][iz]     *= vp[ix][iz];
      vn[ix][iz]     *= vn[ix][iz];
      ONE_P_2_DELTA   =  vn[ix][iz] / vp[ix][iz];
      vh[ix][iz]      =  vn[ix][iz] * (1.0f + 2.0f * eta[ix][iz]);
      lin_eta[ix][iz] = eta[ix][iz] * ONE_P_2_DELTA;
      eps    [ix][iz] = ((1.0f + 2.0f * eta[ix][iz]) * 
                          ONE_P_2_DELTA - 1.0f) * 0.5f;     
    }
  }


  float *kx = sf_floatalloc(nkx);
  float *kz = sf_floatalloc(nkz);
  for (int ikx=0; ikx<nkx; ++ikx) {
    kx[ikx] = okx + ikx * dkx;
    if (ikx >= nkx/2) kx[ikx] = (nkx - ikx) * dkx;
    kx[ikx] *= 2 * SF_PI;
    kx[ikx] *= kx[ikx];
  }
  for (int ikz=0; ikz<nkz; ++ikz) {
    kz[ikz] = okz + ikz * dkz;
    kz[ikz] *= 2 * SF_PI;
    kz[ikz] *= kz[ikz];
  }

  if (adj) {
    ai = ar; ao = as;
  } else {
    ai = as; ao = ar;
  }

  if (opt) {
    sf_oaxa(file_out, ao, 1);
    sf_oaxa(file_out, at, 2);
  } else {
    sf_oaxa(file_out, az, 1);
    sf_oaxa(file_out, ax, 2);
    sf_oaxa(file_out, at, 3);
  }
  sf_fileflush(file_out, NULL);

  if (spt) {
    src2d = pt2dalloc1(ns);
    pt2dread1(file_src, src2d, ns, 2);
    cssinc = sinc2d_make(ns, src2d, nz, nx, dz, dx, oz, ox);
    ws = sf_floatalloc(ns);
    if (adj) { cosinc = cssinc;  wo = ws; }
    else     { cisinc = cssinc;  wi = ws; }
  }
  if (rpt) {
    rec2d = pt2dalloc1(nr);
    pt2dread1(file_rec, rec2d, nr, 2);
    crsinc = sinc2d_make(nr, rec2d, nz, nx, dz, dx, oz, ox);
    wr = sf_floatalloc(nr);
    if (adj) { cisinc = crsinc;  wi = wr; }
    else     { cosinc = crsinc;  wo = wr; }
  }

  u0 = sf_floatalloc2(nz, nx);
  u1 = sf_floatalloc2(nz, nx);
  float *rwave = (float *) fftwf_malloc(nzpad*nxpad*sizeof(float));
  float *rwavem = (float *) fftwf_malloc(nzpad*nxpad*sizeof(float));
  fftwf_complex *cwave = (fftwf_complex *) fftwf_malloc(nkz*nkx*sizeof(fftwf_complex));
  fftwf_complex *cwavem = (fftwf_complex *) fftwf_malloc(nkz*nkx*sizeof(fftwf_complex));
  /* float *rwavem = (float *) fftwf_malloc(nzpad*nxpad*sizeof(float));
  fftwf_complex *cwave = (fftwf_complex *) fftwf_malloc(nkz*nkx*sizeof(fftwf_complex));
  fftwf_complex *cwavem = (fftwf_complex *) fftwf_malloc(nkz*nkx*sizeof(fftwf_complex)); */

  /* boundary conditions */
  float **ucut = NULL;
  float *damp = NULL;
  if (!(ipt &&opt)) ucut = sf_floatalloc2(nz0, nx0);
  damp = damp_make(nb);
    
  float wt = 1./(nxpad * nzpad);
  wt *= dt * dt;
  fftwf_plan forward_plan;
  fftwf_plan inverse_plan;
#ifdef SF_HAS_FFTW_OMP
  fftwf_init_threads();
  fftwf_plan_with_nthreads(omp_get_max_threads());
#endif
  forward_plan = fftwf_plan_dft_r2c_2d(nxpad, nzpad,
              rwave, cwave, FFTW_MEASURE); 
#ifdef SF_HAS_FFTW_OMP
  fftwf_plan_with_nthreads(omp_get_max_threads());
#endif
  inverse_plan = fftwf_plan_dft_c2r_2d(nxpad, nzpad,
              cwavem, rwavem, FFTW_MEASURE); 
  int itb, ite, itc;
  if (adj) {
    itb = nt -1; ite = -1; itc = -1;
  } else {
    itb = 0; ite = nt; itc = 1;
  }

  if (adj) {
    for (int it=0; it<nt; it++) {
      if (opt) sf_floatwrite(wo, sf_n(ao), file_out);
      else     sf_floatwrite(ucut[0], nz0*nx0, file_out);
    }
    sf_seek(file_out, 0, SEEK_SET);
  }

  float **ptrtmp = NULL;
  memset(u0[0], 0, sizeof(float)*nz*nx);
  memset(u1[0], 0, sizeof(float)*nz*nx);
  memset(rwave, 0, sizeof(float)*nzpad*nxpad);
  memset(rwavem, 0, sizeof(float)*nzpad*nxpad);
  memset(cwave, 0, sizeof(float)*nkz*nkx*2);
  memset(cwavem, 0, sizeof(float)*nkz*nkx*2);


  for (int it=itb; it!=ite; it+=itc) { if (verb) sf_warning("it = %d;",it);
    double tic = omp_get_wtime();
    if (ipt) {
      if (adj) sf_seek(file_inp, (off_t)(it)*sizeof(float)*sf_n(ai), SEEK_SET);
      sf_floatread(wi, sf_n(ai), file_inp);
      for (int i=0; i<sf_n(ai); i++)
        wi[i] *= dt* dt;
    } else {
      if (adj) sf_seek(file_inp, (off_t)(it)*sizeof(float)*nz0*nx0, SEEK_SET);
      sf_floatread(ucut[0], nz0*nx0, file_inp);
      for (int j=0; j<nx0; j++)
      for (int i=0; i<nz0; i++)
        ucut[j][i] *= dt * dt;
    }

    /* apply absorbing boundary condition: E \times u@n-1 */
    damp2d_apply(u0, damp, nz, nx, nb);
    fft_stepforward(u0, u1, rwave, rwavem, cwave, cwavem,
        vp, vn, eta, vh, eps, lin_eta, kz, kx,
        forward_plan, inverse_plan,
        nz, nx, nzpad, nxpad, nkz, nkx, wt, adj);

    // sinc2d_inject1(u0, ws[it][s_idx], cssinc[s_idx]);
    if (ipt) sinc2d_inject(u0, wi, cisinc);
    else     wfld2d_inject(u0, ucut, nz0, nx0, nb);

    /* apply absorbing boundary condition: E \times u@n+1 */
    damp2d_apply(u0, damp, nz, nx, nb);

    /* loop over pointers */
    ptrtmp = u0;  u0 = u1;  u1 = ptrtmp;
    
    if (opt) {
      if (adj) sf_seek(file_out, (off_t)(it)*sizeof(float)*sf_n(ao),SEEK_SET);
      sinc2d_extract(u0, wo, cosinc);
      sf_floatwrite(wo, sf_n(ao), file_out);
    } else {
      if (adj) sf_seek(file_out, (off_t)(it)*sizeof(float)*nz0*nx0,SEEK_SET);
      wwin2d(ucut, u0, nz0, nx0, nb);
      sf_floatwrite(ucut[0], nz0*nx0, file_out);
    }

    double toc = omp_get_wtime();
    if (verb) fprintf(stderr," clock = %lf;", toc-tic);
  } /* END OF TIME LOOP */
  return 0;
}

static void
fft_stepforward(
    float **u0, float **u1,
    float *rwave, float *rwavem,
    fftwf_complex *cwave, fftwf_complex *cwavem,
    float **vp, float **vn, float **eta,
    float **vh, float **eps, float **lin_eta,
    float *kz, float *kx,
    fftwf_plan forward_plan, fftwf_plan inverse_plan,
    int nz, int nx, int nzpad, int nxpad,
    int nkz, int nkx,
    float wt, bool adj)

{

#pragma omp parallel for schedule(dynamic,1)
  for (int ix=0; ix<nxpad; ix++) {
    memset(&rwave[ix*nzpad], 0, sizeof(float)*nzpad);
    memset(&rwavem[ix*nzpad], 0, sizeof(float)*nzpad);
    memset(&cwave[ix*nkz], 0, sizeof(fftwf_complex)*nkz);
    memset(&cwavem[ix*nkz], 0, sizeof(fftwf_complex)*nkz);
  }

  if (adj) { /* adjoint modeling */

    /* adj term 1 */ 

#pragma omp parallel for schedule(dynamic,1)
      for (int j=0; j<nx; j++) {
        for (int i=0; i<nz; i++) {
          int jj = j*nzpad + i;
          u0[j][i] = 2.0f * u1[j][i] - u0[j][i]; 
          rwave[jj] = u1[j][i] * vh[j][i];
        }
      }      

    /* --- 2D forward Fourier transform ---*/
    fftwf_execute(forward_plan);

#pragma omp parallel for schedule(dynamic,1)
    for (int ikx=0; ikx<nkx; ++ikx) {
      for (int ikz=0; ikz<nkz; ++ikz) {
        int idx = ikx * nkz + ikz;
        cwavem[idx] += cwave[idx] * kx[ikx];
        //cwavem[idx] =  cwave[idx] * (kx[ikx]+kz[ikz]);
      }
    }

    /* adj term 2 */ 

#pragma omp parallel for schedule(dynamic,1)
      for (int j=0; j<nx; j++) {
        for (int i=0; i<nz; i++) {
          int jj = j*nzpad + i;
          rwave[jj] = u1[j][i] * vp[j][i];
        }
      }

    /* --- 2D forward Fourier transform ---*/
    fftwf_execute(forward_plan);

#pragma omp parallel for schedule(dynamic,1)
    for (int ikx=0; ikx<nkx; ++ikx) {
      for (int ikz=0; ikz<nkz; ++ikz) {
        int idx = ikx * nkz + ikz;
        cwavem[idx] += cwave[idx] * kz[ikz];
      }
    }

    /* adj term 3 */ 

#pragma omp parallel for schedule(dynamic,1)
      for (int j=0; j<nx; j++) {
        for (int i=0; i<nz; i++) {
          int jj = j*nzpad + i;
          rwave[jj] = u1[j][i] * vn[j][i] * eta[j][i] * 2.f;
        }
      }

    /* --- 2D forward Fourier transform ---*/
    fftwf_execute(forward_plan);

#pragma omp parallel for schedule(dynamic,1)
    for (int ikx=0; ikx<nkx; ++ikx) {
      float inv_kx = 1. / kx[ikx];
      for (int ikz=0; ikz<nkz; ++ikz) {
        float inv_kz = 1. / kz[ikz];
        /* float ratio = kx_ * kz_ / (kx_ + kz_); */
        float ratio = 0.f;
        if (isinf(inv_kx) || isinf(inv_kz)) ratio = 0.f;
        else ratio = 1./ (inv_kx + inv_kz);
        /* sf_warning("ratio = %f ", ratio); */
        int idx = ikx * nkz + ikz;
        cwavem[idx] -=  cwave[idx] * ratio;
      }
    }

    /* adj term 4 */ 

#pragma omp parallel for schedule(dynamic,1)
      for (int j=0; j<nx; j++) {
        for (int i=0; i<nz; i++) {
          int jj = j*nzpad + i;
          rwave[jj] = u1[j][i] * vn[j][i] * eps[j][i] * 8.f * Q1;
        }
      }

    /* --- 2D forward Fourier transform ---*/
    fftwf_execute(forward_plan);

#pragma omp parallel for schedule(dynamic,1)
    for (int ikx=0; ikx<nkx; ++ikx) {
      float inv_kx = 1. / kx[ikx];
      for (int ikz=0; ikz<nkz; ++ikz) {
        float inv_kz = 1. / kz[ikz];
        float ratio = 0.f;
        if (isinf(inv_kx) || isinf(inv_kz)) ratio = 0.f;
        else ratio = inv_kz / (inv_kx + inv_kz)*(inv_kx + inv_kz);
        int idx = ikx * nkz + ikz;
        cwavem[idx] += cwave[idx] * ratio;
      }
    }

    /* adj term 5 */ 

#pragma omp parallel for schedule(dynamic,1)
      for (int j=0; j<nx; j++) {
        for (int i=0; i<nz; i++) {
          int jj = j*nzpad + i;
          rwave[jj] = u1[j][i] * vn[j][i] * eta[j][i] * 32.f * Q1 * Q2 *
                                            lin_eta[j][i];
        }
      }

    /* --- 2D forward Fourier transform ---*/
    fftwf_execute(forward_plan);

#pragma omp parallel for schedule(dynamic,1)
    for (int ikx=0; ikx<nkx; ++ikx) {
      float inv_kx = 1. / kx[ikx];
      for (int ikz=0; ikz<nkz; ++ikz) {
        float inv_kz = 1. / kz[ikz];
        float ratio = 0.f;
        if (isinf(inv_kx) || isinf(inv_kz)) ratio = 0.f;
        else ratio = 1./ (inv_kx + inv_kz)*(inv_kx + inv_kz);
        int idx = ikx * nkz + ikz;
        cwavem[idx] -= cwave[idx] * ratio;
      }
    }

    fftwf_execute(inverse_plan);

#pragma omp parallel for schedule(dynamic,1)
    for (int j=0; j<nx; j++) {
      for (int i=0; i<nz; i++) {
        int jj = j*nzpad + i;
        u0[j][i] -= wt * rwavem[jj];
      }
    }

  } else { /* forward modeling */
#pragma omp parallel for schedule(dynamic,1)
      for (int j=0; j<nx; j++) {
        for (int i=0; i<nz; i++) {
          int jj = j*nzpad + i;
          u0[j][i] = 2.0f *u1[j][i] - u0[j][i]; 
          rwave[jj] = u1[j][i];
        }
      }

    /* --- 2D forward Fourier transform ---*/
    fftwf_execute(forward_plan);

    /* term 1 */
#pragma omp parallel for schedule(dynamic,1)
    for (int ikx=0; ikx<nkx; ++ikx) {
      for (int ikz=0; ikz<nkz; ++ikz) {
        int idx = ikx * nkz + ikz;
        cwavem[idx] =  cwave[idx] * kx[ikx];
      }
    }
    fftwf_execute(inverse_plan);
#pragma omp parallel for schedule(dynamic,1)
    for (int j=0; j<nx; j++) {
      for (int i=0; i<nz; i++) {
        int jj = j*nzpad + i;
        u0[j][i] -= wt * rwavem[jj] * vh[j][i];
      }
    }

    /* term 2 */
#pragma omp parallel for schedule(dynamic,1)
    for (int ikx=0; ikx<nkx; ++ikx) {
      for (int ikz=0; ikz<nkz; ++ikz) {
        int idx = ikx * nkz + ikz;
        cwavem[idx] =  cwave[idx] * kz[ikz];
      }
    }
    fftwf_execute(inverse_plan);
#pragma omp parallel for schedule(dynamic,1)
    for (int j=0; j<nx; j++) {
      for (int i=0; i<nz; i++) {
        int jj = j*nzpad + i;
        u0[j][i] -= wt * rwavem[jj] * vp[j][i];
      }
    }
    /* term 3 */
#pragma omp parallel for schedule(dynamic,1)
    for (int ikx=0; ikx<nkx; ++ikx) {
      float inv_kx = 1. / kx[ikx];
      for (int ikz=0; ikz<nkz; ++ikz) {
        float inv_kz = 1. / kz[ikz];
        /* float ratio = kx_ * kz_ / (kx_ + kz_); */
        float ratio = 0.f;
        if (isinf(inv_kx) || isinf(inv_kz)) ratio = 0.f;
        else ratio = 1./ (inv_kx + inv_kz);
        /* sf_warning("ratio = %f ", ratio); */
        int idx = ikx * nkz + ikz;
        cwavem[idx] =  cwave[idx] * ratio;
      }
    }
    fftwf_execute(inverse_plan);
#pragma omp parallel for schedule(dynamic,1)
    for (int j=0; j<nx; j++) {
      for (int i=0; i<nz; i++) {
        int jj = j*nzpad + i;
        u0[j][i] += wt * rwavem[jj] * 2.f * vn[j][i] * eta[j][i];
      }
    }
    /* term 4 */
#pragma omp parallel for schedule(dynamic,1)
    for (int ikx=0; ikx<nkx; ++ikx) {
      float inv_kx = 1. / kx[ikx];
      for (int ikz=0; ikz<nkz; ++ikz) {
        float inv_kz = 1. / kz[ikz];
        float ratio = 0.f;
        if (isinf(inv_kx) || isinf(inv_kz)) ratio = 0.f;
        else ratio = inv_kz / (inv_kx + inv_kz)*(inv_kx + inv_kz);
        int idx = ikx * nkz + ikz;
        cwavem[idx] = cwave[idx] * ratio;
      }
    }
    fftwf_execute(inverse_plan);
#pragma omp parallel for schedule(dynamic,1)
    for (int j=0; j<nx; j++) {
      for (int i=0; i<nz; i++) {
        int jj = j*nzpad + i;
        u0[j][i] -= wt * rwavem[jj] * 8.f * Q1 * 
                    vn[j][i] * eps[j][i];
      }
    }
    /* term 5 */
#pragma omp parallel for schedule(dynamic,1)
    for (int ikx=0; ikx<nkx; ++ikx) {
      float inv_kx = 1. / kx[ikx];
      for (int ikz=0; ikz<nkz; ++ikz) {
        float inv_kz = 1. / kz[ikz];
        float ratio = 0.f;
        if (isinf(inv_kx) || isinf(inv_kz)) ratio = 0.f;
        else ratio = 1./ (inv_kx + inv_kz)*(inv_kx + inv_kz);
        int idx = ikx * nkz + ikz;
        cwavem[idx] = cwave[idx] * ratio;
      }
    }
    fftwf_execute(inverse_plan);
#pragma omp parallel for schedule(dynamic,1)
    for (int j=0; j<nx; j++) {
      for (int i=0; i<nz; i++) {
        int jj = j*nzpad + i;
        u0[j][i] += wt * rwavem[jj] * 32.f * Q1 * Q2 * 
                    vn[j][i] * eta[j][i] * lin_eta[j][i];
      }
    }
  }

  return;
}
