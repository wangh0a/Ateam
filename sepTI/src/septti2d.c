// DONE: TTI case for linear approximation (equation 14), with absorbing boundary
// TODO: adjoint mode
#include <rsf.h>
#include <omp.h>
#include <complex.h>
#include <fftw3.h>
#include "util.h"
#include "sinc.h"
#include <math.h>

#define EPSILON 1e-6
static inline bool is_zero(float a) { return fpclassify(a) == FP_ZERO;}
/* #define is_zero(a) essentiallyEqual(a, 0.f, EPSILON) */

static void
fft_stepforward(
    float **u0, float **u1,
    float *rwave, float *rwavem,
    fftwf_complex *cwave, fftwf_complex *cwavem,
    float **vp, float **eps, float **del, float **st, float **ct, float **st2, float **ct2,
    float *kz, float *kx, float *kz2, float *kx2,
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

  float **vp, **eps, **del, **theta;
  vp  = sf_floatalloc2(nz, nx);
  eps  = sf_floatalloc2(nz, nx);
  del = sf_floatalloc2(nz, nx);
  theta = sf_floatalloc2(nz, nx);
  float **tmparray = sf_floatalloc2(nz0, nx0);
  sf_floatread(tmparray[0], nz0*nx0, file_mdl); expand2d(vp[0], tmparray[0], nz, nx, nz0, nx0);
  sf_floatread(tmparray[0], nz0*nx0, file_mdl); expand2d(eps[0], tmparray[0], nz, nx, nz0, nx0);
  sf_floatread(tmparray[0], nz0*nx0, file_mdl); expand2d(del[0], tmparray[0], nz, nx, nz0, nx0);
  sf_floatread(tmparray[0], nz0*nx0, file_mdl); expand2d(theta[0], tmparray[0], nz, nx, nz0, nx0);

  float **st, **ct, **st2, **ct2;
  st = sf_floatalloc2(nz, nx);
  st2 = sf_floatalloc2(nz, nx);
  ct = sf_floatalloc2(nz, nx);
  ct2 = sf_floatalloc2(nz, nx);
  for (int ix=0; ix<nx; ix++) {
    for (int iz=0; iz<nz; iz++){
      vp[ix][iz] *= vp[ix][iz];
      st[ix][iz] = sinf(theta[ix][iz] * SF_PI / 180.f);
      ct[ix][iz] = cosf(theta[ix][iz] * SF_PI / 180.f);
      st2[ix][iz] = st[ix][iz] * st[ix][iz];
      ct2[ix][iz] = ct[ix][iz] * ct[ix][iz];
      /* sf_warning("theta=%f ct=%f st=%f ct2=%f st2=%f",theta[ix][iz],ct[ix][iz],st[ix][iz],ct2[ix][iz],st2[ix][iz]); */
    }
  }


  float *kx = sf_floatalloc(nkx);
  float *kz = sf_floatalloc(nkz);
  float *kx2 = sf_floatalloc(nkx);
  float *kz2 = sf_floatalloc(nkz);
  for (int ikx=0; ikx<nkx; ++ikx) {
    kx[ikx] = okx + ikx * dkx;
    /* if (ikx >= nkx/2) kx[ikx] = (nkx - ikx) * dkx; */
    if (ikx >= nkx/2) kx[ikx] = (ikx - nkx) * dkx;
    kx[ikx] *= 2 * SF_PI;
    kx2[ikx] = kx[ikx] * kx[ikx];
  }
  for (int ikz=0; ikz<nkz; ++ikz) {
    kz[ikz] = okz + ikz * dkz;
    kz[ikz] *= 2 * SF_PI;
    kz2[ikz] = kz[ikz] * kz[ikz];
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
#ifdef _OPENMP
    double tic = omp_get_wtime();
#endif
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
        vp, eps, del, st, ct, st2, ct2, kz, kx, kz2, kx2,
        /* vp, eps, del, st, ct, st2, ct2, kz, kx, kz2, kx2, */
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
#ifdef _OPENMP
    double toc = omp_get_wtime();
    if (verb) fprintf(stderr," clock = %lf;", toc-tic);
#endif
  } /* END OF TIME LOOP */
  return 0;
}

static void
fft_stepforward(
    float **u0, float **u1,
    float *rwave, float *rwavem,
    fftwf_complex *cwave, fftwf_complex *cwavem,
    float **vp, float **eps, float **del, float **st, float **ct, float **st2, float **ct2,
    float *kz, float *kx, float *kz2, float *kx2,
    fftwf_plan forward_plan, fftwf_plan inverse_plan,
    int nz, int nx, int nzpad, int nxpad,
    int nkz, int nkx,
    float wt, bool adj)

{
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
  for (int ix=0; ix<nxpad; ix++) {
    memset(&rwave[ix*nzpad], 0, sizeof(float)*nzpad);
    memset(&rwavem[ix*nzpad], 0, sizeof(float)*nzpad);
    memset(&cwave[ix*nkz], 0, sizeof(fftwf_complex)*nkz);
    memset(&cwavem[ix*nkz], 0, sizeof(fftwf_complex)*nkz);
  }

  if (adj) { /* adjoint modeling */
    /* add adjoint code here */
  } else { /* forward modeling */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
      for (int j=0; j<nx; j++) {
        for (int i=0; i<nz; i++) {
          int jj = j*nzpad + i;
          u0[j][i] = 2.0f *u1[j][i] - u0[j][i]; 
          rwave[jj] = u1[j][i];
        }
      }

    /* --- 2D forward Fourier transform ---*/
    fftwf_execute(forward_plan);

    /* term 1 ------------- */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
    for (int ikx=0; ikx<nkx; ++ikx) {
      for (int ikz=0; ikz<nkz; ++ikz) {
        int idx = ikx * nkz + ikz;
        cwavem[idx] = cwave[idx] * (kx2[ikx] + kz2[ikz]);
      }
    }
    fftwf_execute(inverse_plan);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
    for (int j=0; j<nx; j++) {
      for (int i=0; i<nz; i++) {
        int jj = j*nzpad + i;
        u0[j][i] -= wt * rwavem[jj] * vp[j][i];
      }
    }

    /* term 2 */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
    for (int ikx=0; ikx<nkx; ++ikx) {
      for (int ikz=0; ikz<nkz; ++ikz) {
        int idx = ikx * nkz + ikz;
        float ratio = 0.f;
        float k2sum = kx2[ikx] + kz2[ikz];
        if (is_zero(k2sum)) ratio = 0.f;
        else ratio  = kx2[ikx] * kx2[ikx] / k2sum;
        cwavem[idx] = cwave[idx] * ratio;
      }
    }
    fftwf_execute(inverse_plan);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
    for (int j=0; j<nx; j++) {
      for (int i=0; i<nz; i++) {
        int jj = j*nzpad + i;
        u0[j][i] -= wt * rwavem[jj] * vp[j][i] 
                    * 2.f * (eps[j][i] * ct2[j][i] * ct2[j][i] + del[j][i] * st2[j][i] * ct2[j][i]); 
      }
    }

    /* term 3 */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
    for (int ikx=0; ikx<nkx; ++ikx) {
      for (int ikz=0; ikz<nkz; ++ikz) {
        int idx = ikx * nkz + ikz;
        float ratio = 0.f;
        float k2sum = kx2[ikx] + kz2[ikz];
        if (!is_zero(k2sum)) ratio  = kz2[ikz] * kz2[ikz] / k2sum;
        cwavem[idx] = cwave[idx] * ratio;
      }
    }
    fftwf_execute(inverse_plan);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
    for (int j=0; j<nx; j++) {
      for (int i=0; i<nz; i++) {
        int jj = j*nzpad + i;
        u0[j][i] -= wt * rwavem[jj] * vp[j][i] 
                   * 2.f * (eps[j][i] * st2[j][i] * st2[j][i] + del[j][i] * st2[j][i] * ct2[j][i]); 
      }
    }

    /* term 4 */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
    for (int ikx=0; ikx<nkx; ++ikx) {
      for (int ikz=0; ikz<nkz; ++ikz) {
        int idx = ikx * nkz + ikz;
        float ratio = 0.f;
        float k2sum = kx2[ikx] + kz2[ikz];
        if (!is_zero(k2sum)) ratio = kx[ikx] * kx2[ikx] * kz[ikz] / k2sum;
        cwavem[idx] = cwave[idx] * ratio;
      }
    }
    fftwf_execute(inverse_plan);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
    for (int j=0; j<nx; j++) {
      for (int i=0; i<nz; i++) {
        int jj = j*nzpad + i;
        u0[j][i] -= wt * rwavem[jj] * vp[j][i] 
                    * (4.f * del[j][i] * st[j][i] * ct[j][i] * (ct2[j][i] - st2[j][i])
                       - 8.f * eps[j][i] * st[j][i] * ct[j][i] * ct2[j][i]);
      }
    }


    /* term 5 */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
    for (int ikx=0; ikx<nkx; ++ikx) {
      for (int ikz=0; ikz<nkz; ++ikz) {
        int idx = ikx * nkz + ikz;
        float ratio = 0.f;
        float k2sum = kx2[ikx] + kz2[ikz];
        if (!is_zero(k2sum)) ratio = kx[ikx] * kz2[ikz] * kz[ikz] / k2sum;
        cwavem[idx] = cwave[idx] * ratio;
      }
    }
    fftwf_execute(inverse_plan);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
    for (int j=0; j<nx; j++) {
      for (int i=0; i<nz; i++) {
        int jj = j*nzpad + i;
        u0[j][i] -= wt * rwavem[jj] * vp[j][i] 
                    * (- 4.f * del[j][i] * st[j][i] * ct[j][i] * (ct2[j][i] - st2[j][i])
                       - 8.f * eps[j][i] * st[j][i] * ct[j][i] * st2[j][i]);
      }
    }


    /* term 6 */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
    for (int ikx=0; ikx<nkx; ++ikx) {
      for (int ikz=0; ikz<nkz; ++ikz) {
        int idx = ikx * nkz + ikz;
        float ratio = 0.f;
        float k2sum = kx2[ikx] + kz2[ikz];
        if (!is_zero(k2sum)) ratio = kx2[ikx] * kz2[ikz] / k2sum;
        cwavem[idx] = cwave[idx] * ratio;
      }
    }
    fftwf_execute(inverse_plan);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
    for (int j=0; j<nx; j++) {
      for (int i=0; i<nz; i++) {
        int jj = j*nzpad + i;
        u0[j][i] -= wt * rwavem[jj] * vp[j][i] 
                    * (2.f * del[j][i] * (ct2[j][i] - st2[j][i]) * (ct2[j][i] - st2[j][i])
                       + 4.f * (3.f * eps[j][i] - del[j][i])* st2[j][i] * ct2[j][i]);
      }
    }

  /* ================================================================== */
#if 1
    /* additional term 1  :   kx^6 */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
    for (int ikx=0; ikx<nkx; ++ikx) {
      for (int ikz=0; ikz<nkz; ++ikz) {
        int idx = ikx * nkz + ikz;
        float ratio = 0.f;
        float k2sum = kx2[ikx] + kz2[ikz];
        k2sum *= k2sum;
        if (!is_zero(k2sum)) ratio = kx2[ikx] * kx2[ikx] * kx2[ikx] / k2sum;
        cwavem[idx] = cwave[idx] * ratio;
      }
    }
    fftwf_execute(inverse_plan);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
    for (int j=0; j<nx; j++) {
      for (int i=0; i<nz; i++) {
        int jj = j*nzpad + i;
        u0[j][i] -= wt * rwavem[jj] * vp[j][i] * 4.f * eps[j][i] * (eps[j][i] - del[j][i])
                    * ct2[j][i] * ct2[j][i] * st2[j][i];
      }
    }
  
    /* additional term 2  :   kz^6 */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
    for (int ikx=0; ikx<nkx; ++ikx) {
      for (int ikz=0; ikz<nkz; ++ikz) {
        int idx = ikx * nkz + ikz;
        float ratio = 0.f;
        float k2sum = kx2[ikx] + kz2[ikz];
        k2sum *= k2sum;
        if (!is_zero(k2sum)) ratio = kz2[ikz] * kz2[ikz] * kz2[ikz] / k2sum;
        cwavem[idx] = cwave[idx] * ratio;
      }
    }
    fftwf_execute(inverse_plan);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
    for (int j=0; j<nx; j++) {
      for (int i=0; i<nz; i++) {
        int jj = j*nzpad + i;
        u0[j][i] -= wt * rwavem[jj] * vp[j][i] * 4.f * eps[j][i] * (eps[j][i] - del[j][i])
                    * ct2[j][i] * st2[j][i] * st2[j][i];
      }
    }

    /* additional term 3  :   kx^5kz */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
    for (int ikx=0; ikx<nkx; ++ikx) {
      for (int ikz=0; ikz<nkz; ++ikz) {
        int idx = ikx * nkz + ikz;
        float ratio = 0.f;
        float k2sum = kx2[ikx] + kz2[ikz];
        k2sum *= k2sum;
        if (!is_zero(k2sum)) ratio = kx2[ikx] * kx2[ikx] * kx[ikx] * kz[ikz] / k2sum;
        cwavem[idx] = cwave[idx] * ratio;
      }
    }
    fftwf_execute(inverse_plan);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
    for (int j=0; j<nx; j++) {
      for (int i=0; i<nz; i++) {
        int jj = j*nzpad + i;
        u0[j][i] -= wt * rwavem[jj] * vp[j][i] * 4.f * eps[j][i] * (eps[j][i] - del[j][i])
                    * 2.f * st[j][i] * ct[j][i] *(ct2[j][i] * ct2[j][i] - 2.f * st2[j][i] * ct2[j][i]);
      }
    }

    /* additional term 4  :   kz^5kx */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
    for (int ikx=0; ikx<nkx; ++ikx) {
      for (int ikz=0; ikz<nkz; ++ikz) {
        int idx = ikx * nkz + ikz;
        float ratio = 0.f;
        float k2sum = kx2[ikx] + kz2[ikz];
        k2sum *= k2sum;
        if (!is_zero(k2sum)) ratio = kz2[ikz] * kz2[ikz] * kz[ikz] * kx[ikx] / k2sum;
        cwavem[idx] = cwave[idx] * ratio;
      }
    }
    fftwf_execute(inverse_plan);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
    for (int j=0; j<nx; j++) {
      for (int i=0; i<nz; i++) {
        int jj = j*nzpad + i;
        u0[j][i] -= wt * rwavem[jj] * vp[j][i] * 4.f * eps[j][i] * (eps[j][i] - del[j][i])
                    * 2.f * st[j][i] * ct[j][i] *(st2[j][i] * st2[j][i] - 2.f * st2[j][i] * ct2[j][i]);
      }
    }


    /* additional term 5  :   kx^4kz^2 */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
    for (int ikx=0; ikx<nkx; ++ikx) {
      for (int ikz=0; ikz<nkz; ++ikz) {
        int idx = ikx * nkz + ikz;
        float ratio = 0.f;
        float k2sum = kx2[ikx] + kz2[ikz];
        k2sum *= k2sum;
        if (!is_zero(k2sum)) ratio = kx2[ikx] * kx2[ikx] * kz2[ikz] / k2sum;
        cwavem[idx] = cwave[idx] * ratio;
      }
    }
    fftwf_execute(inverse_plan);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
    for (int j=0; j<nx; j++) {
      for (int i=0; i<nz; i++) {
        int jj = j*nzpad + i;
        u0[j][i] -= wt * rwavem[jj] * vp[j][i] * 4.f * eps[j][i] * (eps[j][i] - del[j][i])
                    * ct2[j][i] * (ct2[j][i] * ct2[j][i] + 6.f * st2[j][i] * st2[j][i] - 8.f * st2[j][i] * ct2[j][i]);
      }
    }

    /* additional term 6  :   kz^4kx^2 */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
    for (int ikx=0; ikx<nkx; ++ikx) {
      for (int ikz=0; ikz<nkz; ++ikz) {
        int idx = ikx * nkz + ikz;
        float ratio = 0.f;
        float k2sum = kx2[ikx] + kz2[ikz];
        k2sum *= k2sum;
        if (!is_zero(k2sum)) ratio = kz2[ikz] * kz2[ikz] * kx2[ikx] / k2sum;
        cwavem[idx] = cwave[idx] * ratio;
      }
    }
    fftwf_execute(inverse_plan);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
    for (int j=0; j<nx; j++) {
      for (int i=0; i<nz; i++) {
        int jj = j*nzpad + i;
        u0[j][i] -= wt * rwavem[jj] * vp[j][i] * 4.f * eps[j][i] * (eps[j][i] - del[j][i])
                    * st2[j][i] * (st2[j][i] * st2[j][i] + 6.f * ct2[j][i] * ct2[j][i] - 8.f * st2[j][i] * ct2[j][i]);
      }
    }

    /* additional term 7  :   kz^3kx^3 */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
    for (int ikx=0; ikx<nkx; ++ikx) {
      for (int ikz=0; ikz<nkz; ++ikz) {
        int idx = ikx * nkz + ikz;
        float ratio = 0.f;
        float k2sum = kx2[ikx] + kz2[ikz];
        k2sum *= k2sum;
        if (!is_zero(k2sum)) ratio = kz2[ikz] * kz[ikz] * kx2[ikx] * kx[ikx] / k2sum;
        cwavem[idx] = cwave[idx] * ratio;
      }
    }
    fftwf_execute(inverse_plan);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
    for (int j=0; j<nx; j++) {
      for (int i=0; i<nz; i++) {
        int jj = j*nzpad + i;
        u0[j][i] -= wt * rwavem[jj] * vp[j][i] * 4.f * eps[j][i] * (eps[j][i] - del[j][i])
                    * st[j][i] * ct[j][i] * (20.f * st2[j][i] * ct2[j][i] - 4.f);
      }
    }


    /* additional term 8  :   kx^8 + kz^8 */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
    for (int ikx=0; ikx<nkx; ++ikx) {
      float kx8 = kx2[ikx]*kx2[ikx]*kx2[ikx]*kx2[ikx];
      for (int ikz=0; ikz<nkz; ++ikz) {
        int idx = ikx * nkz + ikz;
        float ratio = 0.f;
        float k2sum = kx2[ikx] + kz2[ikz];
        k2sum *= k2sum * k2sum;
        float kz8 = kz2[ikz]*kz2[ikz]*kz2[ikz]*kz2[ikz];
        if (!is_zero(k2sum)) ratio = (kx8 + kz8) / k2sum;
        cwavem[idx] = cwave[idx] * ratio;
      }
    }
    fftwf_execute(inverse_plan);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
    for (int j=0; j<nx; j++) {
      for (int i=0; i<nz; i++) {
        int jj = j*nzpad + i;
        u0[j][i] += wt * rwavem[jj] * vp[j][i] * 4.f * (eps[j][i] - del[j][i]) * (eps[j][i] - del[j][i])
                    * st2[j][i] * st2[j][i] * ct2[j][i] * ct2[j][i];
      }
    }

    /* additional term 9  :   kx^7kz + kz^7kx */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
    for (int ikx=0; ikx<nkx; ++ikx) {
      float kx6 = kx2[ikx]*kx2[ikx]*kx2[ikx];
      for (int ikz=0; ikz<nkz; ++ikz) {
        int idx = ikx * nkz + ikz;
        float ratio = 0.f;
        float k2sum = kx2[ikx] + kz2[ikz];
        k2sum *= k2sum * k2sum;
        float kz6 = kz2[ikz]*kz2[ikz]*kz2[ikz];
        if (!is_zero(k2sum)) ratio = kx[ikx] * kz[ikz] *(kx6 - kz6) / k2sum;
        cwavem[idx] = cwave[idx] * ratio;
      }
    }
    fftwf_execute(inverse_plan);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
    for (int j=0; j<nx; j++) {
      for (int i=0; i<nz; i++) {
        int jj = j*nzpad + i;
        u0[j][i] += wt * rwavem[jj] * vp[j][i] * 4.f * (eps[j][i] - del[j][i]) * (eps[j][i] - del[j][i])
                    * 4.f * st2[j][i] * st[j][i] * ct2[j][i] * ct[j][i] * (st2[j][i] - ct2[j][i]);
      }
    }

    /* additional term 10  :   kx^6kz^2 + kz^6kx^2 */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
    for (int ikx=0; ikx<nkx; ++ikx) {
      float kx4 = kx2[ikx]*kx2[ikx];
      for (int ikz=0; ikz<nkz; ++ikz) {
        int idx = ikx * nkz + ikz;
        float ratio = 0.f;
        float k2sum = kx2[ikx] + kz2[ikz];
        k2sum *= k2sum * k2sum;
        float kz4 = kz2[ikz]*kz2[ikz];
        if (!is_zero(k2sum)) ratio = kx2[ikx] * kz2[ikz] * (kx4 + kz4) / k2sum;
        cwavem[idx] = cwave[idx] * ratio;
      }
    }
    fftwf_execute(inverse_plan);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
    for (int j=0; j<nx; j++) {
      for (int i=0; i<nz; i++) {
        int jj = j*nzpad + i;
        u0[j][i] += wt * rwavem[jj] * vp[j][i] * 4.f * (eps[j][i] - del[j][i]) * (eps[j][i] - del[j][i])
                    * st2[j][i] * ct2[j][i] * (6.f - 28.f * st2[j][i] * ct2[j][i]);
      }
    }

    /* additional term 11  :   kx^5kz^3 + kz^5kx^3 */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
    for (int ikx=0; ikx<nkx; ++ikx) {
      float kx3 = kx2[ikx]*kx[ikx];
      for (int ikz=0; ikz<nkz; ++ikz) {
        int idx = ikx * nkz + ikz;
        float ratio = 0.f;
        float k2sum = kx2[ikx] + kz2[ikz];
        k2sum *= k2sum * k2sum;
        float kz3 = kz2[ikz]*kz[ikz];
        if (!is_zero(k2sum)) ratio = kx3 * kz3 * (kx2[ikx] - kz2[ikz]) / k2sum;
        cwavem[idx] = cwave[idx] * ratio;
      }
    }
    fftwf_execute(inverse_plan);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
    for (int j=0; j<nx; j++) {
      for (int i=0; i<nz; i++) {
        int jj = j*nzpad + i;
        u0[j][i] += wt * rwavem[jj] * vp[j][i] * 4.f * (eps[j][i] - del[j][i]) * (eps[j][i] - del[j][i])
                    * st[j][i] * ct[j][i] * (
                        24.f * st2[j][i] * ct2[j][i] * (st2[j][i] - ct2[j][i])
                        + 4.f * (ct2[j][i]*ct2[j][i]*ct2[j][i] - st2[j][i]*st2[j][i]*st2[j][i])
                      );
      }
    }

    /* additional term 12  :   kx^4kz^4 + kz^4kx^4 */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
    for (int ikx=0; ikx<nkx; ++ikx) {
      for (int ikz=0; ikz<nkz; ++ikz) {
        int idx = ikx * nkz + ikz;
        float ratio = 0.f;
        float k2sum = kx2[ikx] + kz2[ikz];
        k2sum *= k2sum * k2sum;
        if (!is_zero(k2sum)) ratio =  kx2[ikx] * kz2[ikz] * kx2[ikx] * kz2[ikz] / k2sum;
        cwavem[idx] = cwave[idx] * ratio;
      }
    }
    fftwf_execute(inverse_plan);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
    for (int j=0; j<nx; j++) {
      for (int i=0; i<nz; i++) {
        int jj = j*nzpad + i;
        u0[j][i] += wt * rwavem[jj] * vp[j][i] * 4.f * (eps[j][i] - del[j][i]) * (eps[j][i] - del[j][i])
                    * (st2[j][i] * st2[j][i] * (st2[j][i] - 4.f*ct2[j][i]) * (st2[j][i] - 4.f*ct2[j][i])
                     + ct2[j][i] * ct2[j][i] * (ct2[j][i] - 4.f*st2[j][i]) * (ct2[j][i] - 4.f*st2[j][i])  
                     + 4.f * st2[j][i] * ct2[j][i] * (5.f * st2[j][i] * ct2[j][i] - 2.f)
                    );
      }
    }


#endif
  }
  return;
}
