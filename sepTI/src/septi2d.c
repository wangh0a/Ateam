#include <rsf.h>
#include <omp.h>
#include <complex.h>
#include <fftw3.h>
#include "util.h"
#include "sinc.h"

static void 
fft_stepforward(
    float **u0, float **u1,
    float *rwavem, fftwf_complex *cwave, fftwf_complex *cwavem,
    float **lft, float **rht,
    fftwf_plan forward_plan, fftwf_plan inverse_plan,
    int nz, int nx, int nzpad, int nxpad,
    int nkz, int nkx,
    int nrank, float wt, bool adj);

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
  int nzx, nkzkx, nrank;
  
  float **u1, **u0;
  float *ws, *wr;
  float **lft, **rht;
  
  sf_file file_lft = NULL, file_rht = NULL;
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
  file_lft = sf_input("left");
  file_rht = sf_input("right");
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
  float *rwavem = (float *) fftwf_malloc(nzpad*nxpad*sizeof(float));
  fftwf_complex *cwave = (fftwf_complex *) fftwf_malloc(nkz*nkx*sizeof(fftwf_complex));
  fftwf_complex *cwavem = (fftwf_complex *) fftwf_malloc(nkz*nkx*sizeof(fftwf_complex));

  if (!sf_histint(file_lft, "n1", &nzx) || nzx != nz*nx) sf_error("Need n1=%d in left",nz*nx);
  if (!sf_histint(file_lft, "n2", &nrank)) sf_error("Need n2= in left"); 
  if (!sf_histint(file_rht, "n2", &nkzkx) || nkzkx != nkz*nkx) sf_error("Need n2=%d in right",nkz*nkx);

  lft = sf_floatalloc2(nzx, nrank);
  // rht = sf_floatalloc2(nrank, nkzkx);
  rht = sf_floatalloc2(nkzkx, nrank);
  
  sf_floatread(lft[0], nzx*nrank, file_lft);
  sf_floatread(rht[0], nrank*nkzkx, file_rht);
  /* transpose(rht[0], nrank, nkzkx); */

  /* boundary conditions */
  float **ucut = NULL;
  float *damp = NULL;
  if (!(ipt &&opt)) ucut = sf_floatalloc2(nz0, nx0);
  damp = damp_make(nb);
    
  float wt = 1./(nxpad * nzpad);
  fftwf_plan forward_plan;
  fftwf_plan inverse_plan;
#ifdef SF_HAS_FFTW_OMP
  fftwf_init_threads();
  fftwf_plan_with_nthreads(omp_get_max_threads());
#endif
  forward_plan = fftwf_plan_dft_r2c_2d(nxpad, nzpad,
              rwavem, cwave, FFTW_MEASURE); 
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
  memset(cwave, 0, sizeof(float)*nkz*nkx*2);
  memset(cwavem, 0, sizeof(float)*nkz*nkx*2);
  for (int it=itb; it!=ite; it+=itc) { if (verb) sf_warning("it = %d;",it);
    double tic = omp_get_wtime();
    if (ipt) {
      if (adj) sf_seek(file_inp, (off_t)(it)*sizeof(float)*sf_n(ai), SEEK_SET);
      sf_floatread(wi, sf_n(ai), file_inp);
    } else {
      if (adj) sf_seek(file_inp, (off_t)(it)*sizeof(float)*nz0*nx0, SEEK_SET);
      sf_floatread(ucut[0], nz0*nx0, file_inp);
    }

    /* apply absorbing boundary condition: E \times u@n-1 */
    damp2d_apply(u0, damp, nz, nx, nb);
    fft_stepforward(u0, u1, rwavem, cwave, cwavem,
        lft, rht, forward_plan, inverse_plan,
        nz, nx, nzpad, nxpad, nkz, nkx, nrank, wt, adj);

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
    float *rwavem, fftwf_complex *cwave, fftwf_complex *cwavem,
    float **lft, float **rht,
    fftwf_plan forward_plan, fftwf_plan inverse_plan,
    int nz, int nx, int nzpad, int nxpad,
    int nkz, int nkx,
    int nrank, float wt, bool adj)
{
#pragma omp parallel for schedule(dynamic,1)
  for (int ix=0; ix<nxpad; ix++) {
    memset(&rwavem[ix*nzpad], 0, sizeof(float)*nzpad);
    memset(&cwave[ix*nkz], 0, sizeof(fftwf_complex)*nkz);
    memset(&cwavem[ix*nkz], 0, sizeof(fftwf_complex)*nkz);
  }

  if (adj) { /* adjoint modeling */
    for (int im=0; im<nrank; im++) {
#pragma omp parallel for schedule(dynamic,1)
      /* rwavem = L^T_i \schur_dot rwave */
        for (int j=0; j<nx; j++) {
          for (int i=0; i<nz; i++) {
            int ii = j*nz+i;
            int jj = j*nzpad+i;
            rwavem[jj] = lft[im][ii] * u1[j][i];
          }
        }
      /* --- 3D forward Fourier transform ---*/
      fftwf_execute(forward_plan);

#pragma omp parallel for schedule(dynamic,1)
      /* cwavem += R^T_i \schur_dot cwave */
      for (int j=0; j<nkx; j++) {
        for (int ii=0; ii<nkz; ii++) {
          int idx = j * nkz + ii;
          cwavem[idx] += rht[im][idx] * cwave[idx];
        }
      }
    }
    /* --- 3D backward Fourier transform ---*/
    fftwf_execute(inverse_plan);

#pragma omp parallel for schedule(dynamic,1)
      for (int j=0; j<nx; j++) {
        /* FFT centering back for first/second axis */
        for (int i=0; i<nz; i++) {
          int jj = j*nzpad + i;
          u0[j][i] = 2.0f *u1[j][i] - u0[j][i]; 
          /* FFT normalization */
          u0[j][i] += rwavem[jj]*wt;
        }
      }

  } else { /* forward modeling */
#pragma omp parallel for schedule(dynamic,1)
      for (int j=0; j<nx; j++) {
        for (int i=0; i<nz; i++) {
          int jj = j*nzpad + i;
          u0[j][i] = 2.0f *u1[j][i] - u0[j][i]; 
          rwavem[jj] = u1[j][i];
        }
      }

    /* --- 3D forward Fourier transform ---*/
    fftwf_execute(forward_plan);

    for (int im=0; im<nrank; im++) {
      /* element-wise vector multiplication: u@t(kz,kx) * M3(im,:) */
#pragma omp parallel for schedule(dynamic,1)
        for (int j=0; j<nkx; j++) {
          for (int ii=0; ii<nkz; ii++) {
            int idx = j * nkz + ii;
            cwavem[idx] = rht[im][idx] * cwave[idx];
          }
        }

      /* --- 3D backward Fourier transform ---*/
      fftwf_execute(inverse_plan);

    /* element-wise vector multiplication: M1(:,im) * u@t(z,x) */
#pragma omp parallel for schedule(dynamic,1)
        for (int j=0; j<nx; j++) {
          for (int i=0; i<nz; i++) {
            int ii = j*nz + i;
            int jj = j*nzpad + i;
            /* FFT normalization \times wt */
            u0[j][i] += lft[im][ii] * rwavem[jj] * wt;
          }
        }
    }
  }
  return;
}