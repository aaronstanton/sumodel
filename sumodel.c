/* Copyright (c) Aaron Stanton 2013.*/
/* All rights reserved.                       */
/* sufdacoustic  :  $Date: April 2013- Last version Ap    2013  */

#include "su.h"
#include "cwp.h"
#include "segy.h"
#include "header.h"
#include "fftw3.h"
#include <time.h>

#ifndef MARK
#define MARK fprintf(stderr,"%s @ %u\n",__FILE__,__LINE__);fflush(stderr);
#endif

/*********************** self documentation **********************/
char *sdoc[] = {
  " 	   							       ",
  " SUFDACOUSTIC  Acoustic Finite difference modelling                 ",
  "               (2nd order in time, 4th order in space)              ",
  "                                                                    ",
  "           User provides:                                           ",
  "                                                                    ",
  "           Other parameters:                                        ",
  "                                                                    ",
  "           Example coding:                                          ",
  "                                                                    ",
 NULL};
/* Credits:
 * Aaron Stanton
 * Trace header fields accessed: ns, dt, otrav.
 * Last changes: April : 2013 
 */
/**************** end self doc ***********************************/

void fd_step(float **u1, float **u2, float **u3, float **v, float **rho, float **rhoinv, int nz, int nx, float dt, float dz, float dx, float d2z, float d2x, int order);
void pspec_step(float **u1, float **u2, float **u3, float **v, float **rho, float **rhoinv, int nz, int nx, float dt, float dz, float dx);
float fd_approx_deriv1_order1(float f1, float f3, float dx);
float fd_approx_deriv2_order2(float f1, float f2, float f3, float d2x);
float fd_approx_deriv2_order4(float f1, float f2, float f3, float f4, float f5, float d2x);
void u_to_p(float **p, float **u2, float **v, float **rho, float **rhoinv, int nz, int nx, float dt, float dz, float dx);
void ricker_wavelet(float *w, float f,float dt);

int main(int argc, char **argv)
{
  int verbose;
  time_t start,finish;
  double elapsed_time;
  
  int     order, method, nt, nz, nx, it, iz, ix, isz, isx, nsdur;
  int     ih, igz;
  float   dt, dz, dx, d2z, d2x, tmax, zmin, xmin, zmax, xmax, sz, sx, sf, sdur, gz;
  float  rhs;
  float   *w;
  float  **v;
  float **rho;
  float **rhoinv;
  float **u1;
  float **u2;
  float **u3;
  float **dout;
  float ***u;
  cwp_String out;
  segy tr;
  FILE* fpout;   /* file pointer to output file */

  fprintf(stderr,"*******SUFDACOUSTIC*********\n");
  /* Initialize */
  initargs(argc, argv);
  requestdoc(1);
  start=time(0);    

  /* Get parameters */
  if (!getparint("verbose", &verbose))  verbose = 0;
  if (!getparint("method", &method))  method = 1; /* 1 is FD method, otherwise PSUEDOSPECTRAL method is used.*/
  if (!getparint("order", &order))  order = 4;
  if (!getparstring("out",&out)) err("out required."); 
  if (!getparfloat("dt", &dt))  dt = 0.001;
  if (!getparfloat("dz", &dz))  dz = 5;
  if (!getparfloat("dx", &dx))  dx = 5;
  if (!getparfloat("tmax", &tmax))  tmax = 1;
  if (!getparfloat("zmin", &zmin))  zmin = 0;
  if (!getparfloat("xmin", &xmin))  xmin = 0;
  if (!getparfloat("zmax", &zmax))  zmax = 1000;
  if (!getparfloat("xmax", &xmax))  xmax = 1000;
  if (!getparfloat("sz", &sz))  sz = 500;
  if (!getparfloat("sx", &sx))  sx = 500;
  if (!getparfloat("sf", &sf))  sf = 20;
  if (!getparfloat("sdur", &sdur))  sdur = 0.5;
  if (!getparfloat("gz", &gz))  gz = 10;

  d2z = dz*dz;
  d2x = dx*dx;
    
  nt = (int) tmax/dt + 1;
  nz = (int) (zmax-zmin)/dz + 1;
  nx = (int) (xmax-xmin)/dx + 1;
  isz = (int) (sz - zmin)/dz; 
  isx = (int) (sx - xmin)/dx;
  nsdur = (int) trunc(sdur/dt) + 1;
  igz = (int) (gz - zmin)/dz;
 
 
  w  = ealloc1float(nsdur);
  for (it=0;it<nsdur;it++){
    w[it] = 0;
  }  
  ricker_wavelet(w, sf, dt);

  v  = ealloc2float(nz,nx);
  rho= ealloc2float(nz,nx);
  rhoinv= ealloc2float(nz,nx);
  u1 = ealloc2float(nz,nx);
  u2 = ealloc2float(nz,nx);
  u3 = ealloc2float(nz,nx);
  dout = ealloc2float(nt,nx);
  u  = ealloc3float(nt,nz,nx);

  for (iz=0;iz<nz;iz++){  
    for (ix=0;ix<nx;ix++){
      v[ix][iz]=1500; /* setting velocity constant to 1500m/s */  
      rho[ix][iz]=1000; /* setting density constant to 1000kg/m^3 */
      rhoinv[ix][iz]=1/rho[ix][iz]; 
      u1[ix][iz]=0;
      u2[ix][iz]=0;
      u3[ix][iz]=0;
      for (it=0;it<nt;it++){  
      	u[ix][iz][it]=0;
      }
    }
  }
  
  it = 0;
  while (it < nt){
  		
	if (method==1){ 
	  fd_step(u1,u2,u3,v,rho,rhoinv,nz,nx,dt,dz,dx,d2z,d2x,order);
	}
	else{
	  pspec_step(u1,u2,u3,v,rho,rhoinv,nz,nx,dt,dz,dx);
	}
    
    /* all boundaries fixed */
    /* set edges to zero */  
    for (iz=0;iz<nz;iz++){  
      u3[1][iz]    = 0; 
      u3[nx-1][iz] = 0;
    }
    for (ix=0;ix<nx;ix++){  
      u3[ix][1]    = 0; 
      u3[ix][nz-1] = 0;
    }
   
    if (it<nsdur){
      u3[isx][isz] = w[it];
    }  
    for (iz=0;iz<nz;iz++){  
      for (ix=0;ix<nx;ix++){  
        u1[ix][iz]=u2[ix][iz];
        u2[ix][iz]=u3[ix][iz];
        u[ix][iz][it]=u2[ix][iz];
      }
    }
    for (ix=0;ix<nx;ix++){  
      dout[ix][it] = u[ix][igz][it];
    }
    it = it+1;    
  }    
 
  free1float(w);
  free2float(v);
  free2float(rho);
  free2float(rhoinv);
  free2float(u1);
  free2float(u2);
  free2float(u3);
  free3float(u);

  /* ***********************************************************************
  outputting data
  *********************************************************************** */
  fpout = efopen(out,"w");
  for (ih=0;ih<nx;ih++){ 
    /* fprintf(stderr,"ih=%d\n",ih); */
    memcpy((void *) tr.data,(const void *) dout[ih],nt*sizeof(float));
    tr.sx  = (int) sx;
    tr.gx  = (int) xmin + dx*ih;
    tr.ns = nt;
    tr.dt = NINT(dt*1000000.);
    tr.tracl = tr.tracr = ih + 1;
    fputtr(fpout,&tr);
  }
  fclose(fpout);
  /* ***********************************************************************
  end outputting data
  *********************************************************************** */
  free2float(dout);
  
  /******** end of program **********/
  finish=time(0);
  elapsed_time=difftime(finish,start);
  fprintf(stderr,"Total time required: %6.2f \n", elapsed_time);
  
  return EXIT_SUCCESS;
}

void fd_step(float **u1, float **u2, float **u3, float **v, float **rho, float **rhoinv, int nz, int nx, float dt, float dz, float dx, float d2z, float d2x, int order)
{
  int iz, ix;
  float rhs;
  
  if (order==2){
    for (iz=1;iz<nz-1;iz++){ 
      for (ix=1;ix<nx-1;ix++){
		rhs = (rho[ix][iz]*v[ix][iz]*v[ix][iz])*(fd_approx_deriv1_order1(rhoinv[ix][iz-1],rhoinv[ix][iz+1],dz)*fd_approx_deriv1_order1(u2[ix][iz-1],u2[ix][iz+1],dz) +
						 fd_approx_deriv1_order1(rhoinv[ix-1][iz],rhoinv[ix+1][iz],dx)*fd_approx_deriv1_order1(u2[ix-1][iz],u2[ix+1][iz],dx) + 
						 (rhoinv[ix][iz])*fd_approx_deriv2_order2(u2[ix][iz-1],u2[ix][iz],u2[ix][iz+1],d2z) +
						 (rhoinv[ix][iz])*fd_approx_deriv2_order2(u2[ix-1][iz],u2[ix][iz],u2[ix+1][iz],d2x));
		u3[ix][iz] = (dt*dt)*rhs + 2*u2[ix][iz] - u1[ix][iz];
      }
    }
  }
  else { 
    for (iz=2;iz<nz-2;iz++){  
      for (ix=2;ix<nx-2;ix++){  
		rhs = (rho[ix][iz]*v[ix][iz]*v[ix][iz])*(fd_approx_deriv1_order1(rhoinv[ix][iz-1],rhoinv[ix][iz+1],dz)*fd_approx_deriv1_order1(u2[ix][iz-1],u2[ix][iz+1],dz) +
						 fd_approx_deriv1_order1(rhoinv[ix-1][iz],rhoinv[ix+1][iz],dx)*fd_approx_deriv1_order1(u2[ix-1][iz],u2[ix+1][iz],dx) +
						 (rhoinv[ix][iz])*fd_approx_deriv2_order4(u2[ix][iz-2],u2[ix][iz-1],u2[ix][iz],u2[ix][iz+1],u2[ix][iz+2],d2z) +  
						 (rhoinv[ix][iz])*fd_approx_deriv2_order4(u2[ix-2][iz],u2[ix-1][iz],u2[ix][iz],u2[ix+1][iz],u2[ix+2][iz],d2x));
		u3[ix][iz] = (dt*dt)*rhs + 2*u2[ix][iz] - u1[ix][iz];
      }
    }  
  }
    
  return;
}

void pspec_step(float **u1, float **u2, float **u3, float **v, float **rho, float **rhoinv, int nz, int nx, float dt, float dz, float dx)
{

  int nkz, nkx, ikz, ikx, iz, ix;
  float **rhs1;
  float *kz;
  float *kx;
  float *in1a;
  float *in1b;
  complex *out1a;
  complex *out1b;
  complex *in2a;
  complex *in2b;
  float *out2a;
  float *out2b;
  fftwf_plan pfwd_a;
  fftwf_plan prv_a;
  fftwf_plan pfwd_b;
  fftwf_plan prv_b;

  nkz = (int) trunc(nz/2) + 1;
  nkx = (int) trunc(nx/2) + 1;

  rhs1  = ealloc2float(nz,nx);
  kz    = ealloc1float(nkz);
  kx    = ealloc1float(nkx);
  in1a  = ealloc1float(nz);
  in1b  = ealloc1float(nx);
  out1a = ealloc1complex(nkz);
  out1b = ealloc1complex(nkx);
  in2a  = ealloc1complex(nkz);
  in2b  = ealloc1complex(nkx);
  out2a = ealloc1float(nz);
  out2b = ealloc1float(nx);

  for (ikz=0;ikz<nkz;ikz++) kz[ikz]= PI*ikz/(nkz)/dz;
  for (ikx=0;ikx<nkx;ikx++) kx[ikx]= PI*ikx/(nkx)/dx;

/* START routine to approximate second derivative along the z direction */
  pfwd_a = fftwf_plan_dft_r2c_1d(nz, in1a, (fftwf_complex*)out1a, FFTW_ESTIMATE);
  prv_a = fftwf_plan_dft_c2r_1d(nz, (fftwf_complex*)in2a, out2a, FFTW_ESTIMATE);
  for (ix=0;ix<nx;ix++){
    for(iz=0;iz<nz;iz++){
      in1a[iz] = u2[ix][iz];
	}
    fftwf_execute(pfwd_a);  
    for(ikz=0;ikz<nkz;ikz++){
      in2a[ikz].r = -out1a[ikz].i*kz[ikz]; 
      in2a[ikz].i =  out1a[ikz].r*kz[ikz]; 
    }
    fftwf_execute(prv_a); 
    for(iz=0;iz<nz;iz++){
      in1a[iz] = rhoinv[ix][iz]*out2a[iz]/nz; 
	}
    fftwf_execute(pfwd_a);  
    for(ikz=0;ikz<nkz;ikz++){
      in2a[ikz].r = -out1a[ikz].i*kz[ikz]; 
      in2a[ikz].i =  out1a[ikz].r*kz[ikz]; 
    }
    fftwf_execute(prv_a); 
    for(iz=0;iz<nz;iz++){
      rhs1[ix][iz] = out2a[iz]/nz;  
    }
  }
/* END routine to approximate second derivative along the z direction */

/* START routine to approximate second derivative along the x direction */
  pfwd_b = fftwf_plan_dft_r2c_1d(nx, in1b, (fftwf_complex*)out1b, FFTW_ESTIMATE);
  prv_b  = fftwf_plan_dft_c2r_1d(nx, (fftwf_complex*)in2b, out2b, FFTW_ESTIMATE);
  for (iz=0;iz<nz;iz++){
    for(ix=0;ix<nx;ix++){
      in1b[ix] = u2[ix][iz];
	}
    fftwf_execute(pfwd_b);  
    for(ikx=0;ikx<nkx;ikx++){
      in2b[ikx].r = -out1b[ikx].i*kx[ikx];
      in2b[ikx].i =  out1b[ikx].r*kx[ikx]; 
    }
    fftwf_execute(prv_b); 
    for(ix=0;ix<nx;ix++){
      in1b[ix] = rhoinv[iz][ix]*out2b[ix]/nx;
	}
    fftwf_execute(pfwd_b);  
    for(ikx=0;ikx<nkx;ikx++){
      in2b[ikx].r = -out1b[ikx].i*kx[ikx]; 
      in2b[ikx].i =  out1b[ikx].r*kx[ikx]; 
    }
    fftwf_execute(prv_b); 
    for(ix=0;ix<nx;ix++){ 
      u3[ix][iz] = (dt*dt)*((rho[ix][iz]*v[ix][iz]*v[ix][iz])*(rhs1[ix][iz] + out2b[ix]/nx)) + 2*u2[ix][iz] - u1[ix][iz]; 
    }
  }
/* END routine to approximate second derivative along the x direction */

  free2float(rhs1);
  free1float(kz);
  free1float(kx);
  fftwf_destroy_plan(pfwd_a);
  fftwf_destroy_plan(prv_a);
  fftwf_destroy_plan(pfwd_b);
  fftwf_destroy_plan(prv_b);
  fftwf_free(in1a); fftwf_free(out1a);
  fftwf_free(in1b); fftwf_free(out1b);
  fftwf_free(in2a); fftwf_free(out2a);
  fftwf_free(in2b); fftwf_free(out2b);

  return;
}

void u_to_p(float **p, float **u2, float **v, float **rho, int nz, int nx, float dz, float dx)
{
  int iz, ix;
  for (iz=1;iz<nz-1;iz++){ 
    for (ix=1;ix<nx-1;ix++){
	  p[ix][iz] = (-rho[ix][iz]*v[ix][iz]*v[ix][iz])*(fd_approx_deriv1_order1(u2[ix][iz-1],u2[ix][iz+1],dz) + 
	  								                  fd_approx_deriv1_order1(u2[ix-1][iz],u2[ix+1][iz],dx));
    }
  }
  for (iz=0;iz<nz;iz++){ 
	  p[0][iz]  = 0;
	  p[nx][iz] = 0;
  }
  for (ix=0;ix<nx;ix++){ 
	  p[ix][0]  = 0;
	  p[ix][nz] = 0;
  }
  return;
}

float fd_approx_deriv1_order1(float f1, float f3, float dx)
{
  /* 1st order finite difference approximation to 1st derivative*/
  /* see www.en.wikipedia.org/wiki/Finite_difference_coefficients for description*/
  float dfdx;
  dfdx = (f3 - f1)/(2*dx);
  return dfdx;
}

float fd_approx_deriv2_order2(float f1, float f2, float f3, float d2x)
{
  /* 2nd order finite difference approximation to 2nd derivative*/
  /* see www.en.wikipedia.org/wiki/Finite_difference_coefficients for description*/
  float d2fd2x;
  d2fd2x = (f3 - 2*f2 + f1)/d2x;
  return d2fd2x;
}

float fd_approx_deriv2_order4(float f1, float f2, float f3, float f4, float f5, float d2x)
{
  /* 4th order finite difference approximation to 2nd derivative*/
  /* see www.en.wikipedia.org/wiki/Finite_difference_coefficients for description*/
  float d2fd2x;
  d2fd2x = (-f5 + 16*f4 - 30*f3 + 16*f2 - f1)/(12*d2x);
  return d2fd2x;
}

void ricker_wavelet(float *w, float f,float dt)
{
  int iw, nw, nc;
  float alpha, beta; 
    
  nw = (int) 2*trunc((float) (2.2/f/dt)/2) + 1;
  nc = (int) trunc((float) nw/2);
 
  for (iw=0;iw<nw-2;iw++){
    alpha = (nc-iw+1)*f*dt*PI;
  	beta = alpha*alpha;
    w[iw] = (1-beta*2)*exp(-beta);
  }

  return;
}
