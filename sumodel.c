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
  " SUMODEL   Modelling seismic data                ",
  "               (FD (method=1) can be 2nd order in time, 2nd order in space, order=2)              ",
  "               (FD (method=1) can be 2nd order in time, 4th order in space, order=4)              ",
  "               (Pseudospectral  is activated using method=2                                       ",
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
 * Trace header fields accessed: 
 * Last changes: April : 2013 
 */
/**************** end self doc ***********************************/

#define	C1	1.125
#define	C2	-0.04166667


void fd_step_acoustic(float **u1, float **u2, float **u3, float **v, float **rho, float **rhoinv, int nz, int nx, float dt, float dz, float dx, float d2z, float d2x, int order);
void  fd_step_elastic_iso(float **ux1, float **uy1, float **uz1, float **ux2, float **uy2, float **uz2, float **ux3, float **uy3, float **uz3, float **v, float **rho, float **rhoinv, int nz, int nx, float dt, float dz, float dx, float d2z, float d2x, int order);
void  pspec_step(float **u1, float **u2, float **u3, float **v, float **rho, float **rhoinv, int nz, int nx, float dt, float dz, float dx);

void update_velocity(int nz, int nx, 
                float **vx, float **vz, float **txx, float **tzz, float **txz,
	            float **rhoinv, float dtx, float dtz);

void update_stress_iso(int nz, int nx, 
	                   float **vx, float **vz, float **txx, float **tzz, float **txz,
	                   float **c11, float **c55, float dtx, float dtz);
	                   	
float fd_approx_deriv1_order2(float f1, float f3, float dx);
float fd_approx_deriv2_order2(float f1, float f2, float f3, float d2x);
float fd_approx_deriv2_order4(float f1, float f2, float f3, float f4, float f5, float d2x);
void v_to_p(float **p, float **vx, float **vy, float **vz, int nz, int nx, float dt, float dz, float dx, float vp_water, float rho_water);
void  ricker_wavelet(float *w, float f,float dt);

int main(int argc, char **argv)
{
  int verbose;
  time_t start,finish;
  double elapsed_time;
  
  int     outcomp, order, method, nt, nz, nx, it, iz, ix, isz, isx, nsdur;
  int     ih, igz;
  float   dt, dz, dx, d2z, d2x, tmax, zmin, xmin, zmax, xmax, sz, sx, sf, sdur, gz;
  float  vp_water,rho_water,dtx,dtz;
  float   *w;
  float  **vp;
  float  **vs;
  float **rho;
  float **c11;
  float **c55;
  float **txx; 
  float **tzz; 
  float **txz;
  float **rhoinv;
  float **vx1;
  float **vy1;
  float **vz1;
  float **p;
  float **dxout;
  float **dyout;
  float **dzout;
  float **dpout;
  float ***vx;
  float ***vy;
  float ***vz;
  cwp_String outz;
  cwp_String outx;
  cwp_String outp;
  segy tr;
  FILE* fpoutz;   /* file pointer to output z file */
  FILE* fpoutx;   /* file pointer to output x file */
  FILE* fpoutp;   /* file pointer to output p file */

  fprintf(stderr,"*******SUMODEL*********\n");
  /* Initialize */
  initargs(argc, argv);
  requestdoc(1);
  start=time(0);    

  /* Get parameters */
  if (!getparint("verbose", &verbose))  verbose = 0;
  if (!getparint("method", &method))  method = 1; /* 1 is FD method, otherwise PSUEDOSPECTRAL method is used.*/
  if (!getparint("order", &order))  order = 4;
  if (!getparint("outcomp", &outcomp))  outcomp = 0;
  if (!getparstring("outz",&outz)) err("outz required."); 
  if (!getparstring("outx",&outx)) err("outx required."); 
  if (!getparstring("outp",&outp)) err("outp required."); 
  if (!getparfloat("dt", &dt))  dt = 0.001;
  if (!getparfloat("dz", &dz))  dz = 5;
  if (!getparfloat("dx", &dx))  dx = 5;
  if (!getparfloat("tmax", &tmax))  tmax = 1;
  if (!getparfloat("zmin", &zmin))  zmin = 0;
  if (!getparfloat("xmin", &xmin))  xmin = 0;
  if (!getparfloat("zmax", &zmax))  zmax = 1000;
  if (!getparfloat("xmax", &xmax))  xmax = 1000;
  if (!getparfloat("sz", &sz))  sz = 200;
  if (!getparfloat("sx", &sx))  sx = 500;
  if (!getparfloat("sf", &sf))  sf = 20;
  if (!getparfloat("sdur", &sdur))  sdur = 0.5;
  if (!getparfloat("gz", &gz))  gz = 100;
  if (!getparfloat("vp_water", &vp_water))  vp_water = 1500;
  if (!getparfloat("rho_water", &rho_water))  rho_water = 1000;

  d2z = dz*dz;
  d2x = dx*dx;
    
  nt = (int) tmax/dt + 1;
  nz = (int) (zmax-zmin)/dz + 1;
  nx = (int) (xmax-xmin)/dx + 1;
  isz = (int) (sz - zmin)/dz; 
  isx = (int) (sx - xmin)/dx;
  nsdur = (int) trunc(sdur/dt) + 1;
  igz = (int) (gz - zmin)/dz;
  dtx = dt/dx; 
  dtz = dt/dz;
 
  w  = ealloc1float(nsdur);
  for (it=0;it<nsdur;it++){
    w[it] = 0;
  }  
  ricker_wavelet(w, sf, dt);

  vp  = ealloc2float(nz,nx);
  vs  = ealloc2float(nz,nx);
  rho= ealloc2float(nz,nx);
  c11 = ealloc2float(nz,nx);
  c55 = ealloc2float(nz,nx);
  txx = ealloc2float(nz,nx);
  tzz = ealloc2float(nz,nx);
  txz = ealloc2float(nz,nx);
  rhoinv= ealloc2float(nz,nx);
  vx1 = ealloc2float(nz,nx);
  vy1 = ealloc2float(nz,nx);
  vz1 = ealloc2float(nz,nx);
  dxout = ealloc2float(nt,nx);
  dyout = ealloc2float(nt,nx);
  dzout = ealloc2float(nt,nx);
  dpout = ealloc2float(nt,nx);
  vx  = ealloc3float(nt,nz,nx);
  vy  = ealloc3float(nt,nz,nx);
  vz  = ealloc3float(nt,nz,nx);
  p = ealloc2float(nz,nx);

  for (iz=0;iz<nz;iz++){  
    for (ix=0;ix<nx;ix++){
      vp[ix][iz]=1500; /* setting P velocity */  
      vs[ix][iz]=0; /* setting S velocity */  
      rho[ix][iz]=1000; /* setting density constant to 1000kg/m^3 */
      c11[ix][iz]=(rho[ix][iz]*vp[ix][iz]*vp[ix][iz] - rho[ix][iz]*vs[ix][iz]*vs[ix][iz]) + 2*rho[ix][iz]*vs[ix][iz]*vs[ix][iz]; /*  lambda + 2*mu  , where lambda = rho*vp^2- mu*/
      c55[ix][iz]=rho[ix][iz]*vs[ix][iz]*vs[ix][iz]; /*  mu  , where mu = rho*vs^2*/
      txx[ix][iz]=0; 
      tzz[ix][iz]=0; 
      txz[ix][iz]=0; 
      rhoinv[ix][iz]=1/rho[ix][iz]; 
      vx1[ix][iz]=0;
      vy1[ix][iz]=0;
      vz1[ix][iz]=0;
      for (it=0;it<nt;it++){  
      	vx[ix][iz][it]=0;
      	vy[ix][iz][it]=0;
      	vz[ix][iz][it]=0;
      }
    }
  }
  it = 0;
  fprintf(stderr,"\n");	
  while (it < nt){
  	fprintf(stderr,"\r%3.0f%% Complete.",(float) (it+1)/nt * 100);	
  	
  	if (it<nsdur){
      vz1[isx][isz] = w[it];
    } 

	/* if (method==1){  */
	  /*
	  fd_step(ux1,uy1,uz1,ux2,uy2,uz2,ux3,uy3,uz3,v,rho,rhoinv,nz,nx,dt,dz,dx,d2z,d2x,order);
	  fd_step_elastic_iso(ux1,uy1,uz1,ux2,uy2,uz2,ux3,uy3,uz3,v,rho,rhoinv,nz,nx,dt,dz,dx,d2z,d2x,order);
      */
      update_velocity(nz,nx,vx1,vz1,txx,tzz,txz,rhoinv,dtx,dtz);

	/* } */
	/* else{ */
	/*
	  pspec_step(ux1,uy1,uz1,ux2,uy2,uz2,ux3,uy3,uz3,v,rho,rhoinv,nz,nx,dt,dz,dx);
	*/  
	/* } */
    /* set velocity boundary conditions */
    for (iz=0;iz<nz;iz++){  
      vx1[1][iz]    = 0; 
      vy1[1][iz]    = 0; 
      vz1[1][iz]    = 0; 
      vx1[nx-1][iz] = 0;
      vy1[nx-1][iz] = 0;
      vz1[nx-1][iz] = 0;
    }
    for (ix=0;ix<nx;ix++){  
      vx1[ix][1]    = 0; 
      vy1[ix][1]    = 0; 
      vz1[ix][1]    = 0; 
      vx1[ix][nz-1] = 0;
      vy1[ix][nz-1] = 0;
      vz1[ix][nz-1] = 0;
    }
   
    for (iz=0;iz<nz;iz++){  
      for (ix=0;ix<nx;ix++){  
        vx[ix][iz][it]=vx1[ix][iz];
        vy[ix][iz][it]=vy1[ix][iz];
        vz[ix][iz][it]=vz1[ix][iz];
      }
    }
    
    update_stress_iso(nz,nx,vx1,vz1,txx,tzz,txz,c11,c55,dtx,dtz);

    /* if (outcomp>0){ */
      v_to_p(p,vx1,vy1,vz1,nz,nx,dt,dz,dx,vp_water,rho_water);
    /* } */

    for (ix=0;ix<nx;ix++){
      /* if (outcomp){ */
      	dpout[ix][it] = p[ix][igz];
      /* } */
     /* else { */
        dxout[ix][it] = vx[ix][igz][it];	
        dyout[ix][it] = vy[ix][igz][it];	
        dzout[ix][it] = vz[ix][igz][it];	
     /* } */
    }
    it = it+1;    
  }    
  fprintf(stderr,"\n");	
 
  free1float(w);
  free2float(vp);
  free2float(vs);
  free2float(rho);
  free2float(c11);
  free2float(c55);
  free2float(txx); 
  free2float(tzz); 
  free2float(txz);
  free2float(rhoinv);
  free2float(vx1);
  free2float(vy1);
  free2float(vz1);
  free2float(p);
  free3float(vx);
  free3float(vy);
  free3float(vz);

  /* ***********************************************************************
  outputting z component
  *********************************************************************** */
  fpoutz = efopen(outz,"w");
  for (ih=0;ih<nx;ih++){ 
    /* fprintf(stderr,"ih=%d\n",ih); */
    memcpy((void *) tr.data,(const void *) dzout[ih],nt*sizeof(float));
    tr.sx  = (int) sx;
    tr.gx  = (int) xmin + dx*ih;
    tr.ns = nt;
    tr.dt = NINT(dt*1000000.);
    tr.tracl = tr.tracr = ih + 1;
    fputtr(fpoutz,&tr);
  }
  fclose(fpoutz);
  /* ***********************************************************************
  end outputting z component
  *********************************************************************** */
 
  /* ***********************************************************************
  outputting x component
  *********************************************************************** */
  fpoutx = efopen(outx,"w");
  for (ih=0;ih<nx;ih++){ 
    /* fprintf(stderr,"ih=%d\n",ih); */
    memcpy((void *) tr.data,(const void *) dxout[ih],nt*sizeof(float));
    tr.sx  = (int) sx;
    tr.gx  = (int) xmin + dx*ih;
    tr.ns = nt;
    tr.dt = NINT(dt*1000000.);
    tr.tracl = tr.tracr = ih + 1;
    fputtr(fpoutx,&tr);
  }
  fclose(fpoutx);
  /* ***********************************************************************
  end outputting x component
  *********************************************************************** */

  /* ***********************************************************************
  outputting p component
  *********************************************************************** */
  fpoutp = efopen(outp,"w");
  for (ih=0;ih<nx;ih++){ 
    /* fprintf(stderr,"ih=%d\n",ih); */
    memcpy((void *) tr.data,(const void *) dpout[ih],nt*sizeof(float));
    tr.sx  = (int) sx;
    tr.gx  = (int) xmin + dx*ih;
    tr.ns = nt;
    tr.dt = NINT(dt*1000000.);
    tr.tracl = tr.tracr = ih + 1;
    fputtr(fpoutp,&tr);
  }
  fclose(fpoutp);
  /* ***********************************************************************
  end outputting p component
  *********************************************************************** */

  free2float(dxout);
  free2float(dyout);
  free2float(dzout);
  free2float(dpout);
  
  /******** end of program **********/
  finish=time(0);
  elapsed_time=difftime(finish,start);
  fprintf(stderr,"Total time required: %6.2f \n", elapsed_time);
  
  return EXIT_SUCCESS;
}

void fd_step_acoustic(float **u1, float **u2, float **u3, float **v, float **rho, float **rhoinv, int nz, int nx, float dt, float dz, float dx, float d2z, float d2x, int order)
{
  int iz, ix;
  float rhs;
  
  if (order==2){
    for (iz=1;iz<nz-1;iz++){ 
      for (ix=1;ix<nx-1;ix++){
		rhs = (rho[ix][iz]*v[ix][iz]*v[ix][iz])*(fd_approx_deriv1_order2(rhoinv[ix][iz-1],rhoinv[ix][iz+1],dz)*fd_approx_deriv1_order2(u2[ix][iz-1],u2[ix][iz+1],dz) +
						 fd_approx_deriv1_order2(rhoinv[ix-1][iz],rhoinv[ix+1][iz],dx)*fd_approx_deriv1_order2(u2[ix-1][iz],u2[ix+1][iz],dx) + 
						 (rhoinv[ix][iz])*fd_approx_deriv2_order2(u2[ix][iz-1],u2[ix][iz],u2[ix][iz+1],d2z) +
						 (rhoinv[ix][iz])*fd_approx_deriv2_order2(u2[ix-1][iz],u2[ix][iz],u2[ix+1][iz],d2x));
		u3[ix][iz] = (dt*dt)*rhs + 2*u2[ix][iz] - u1[ix][iz];
      }
    }
  }
  else { 
    for (iz=2;iz<nz-2;iz++){  
      for (ix=2;ix<nx-2;ix++){  
		rhs = (rho[ix][iz]*v[ix][iz]*v[ix][iz])*(fd_approx_deriv1_order2(rhoinv[ix][iz-1],rhoinv[ix][iz+1],dz)*fd_approx_deriv1_order2(u2[ix][iz-1],u2[ix][iz+1],dz) +
						 fd_approx_deriv1_order2(rhoinv[ix-1][iz],rhoinv[ix+1][iz],dx)*fd_approx_deriv1_order2(u2[ix-1][iz],u2[ix+1][iz],dx) +
						 (rhoinv[ix][iz])*fd_approx_deriv2_order4(u2[ix][iz-2],u2[ix][iz-1],u2[ix][iz],u2[ix][iz+1],u2[ix][iz+2],d2z) +  
						 (rhoinv[ix][iz])*fd_approx_deriv2_order4(u2[ix-2][iz],u2[ix-1][iz],u2[ix][iz],u2[ix+1][iz],u2[ix+2][iz],d2x));
		u3[ix][iz] = (dt*dt)*rhs + 2*u2[ix][iz] - u1[ix][iz];
      }
    }  
  }
    
  return;
}

void update_velocity(int nz, int nx, 
                float **vx, float **vz, float **txx, float **tzz, float **txz,
	            float **rhoinv, float dtx, float dtz)
{
  int iz,ix;
  float dtxx,dtzz,dtxz,dtzx;
  for (iz=2;iz<nz-2;iz++) {
    for (ix=2;ix<nx-2;ix++) {
	  dtxx = C1*(txx[ix][iz+1]-txx[ix][iz]) + C2*(txx[ix][iz+2]-txx[ix][iz-1]);
	  dtxz = C1*(txz[ix+1][iz]-txz[ix][iz]) + C2*(txz[ix+2][iz]-txz[ix-1][iz]);
	  vx[ix][iz] = vx[ix][iz] + ( dtx*dtxx + dtz*dtxz ) * rhoinv[ix][iz];
	  dtzz = C1*(tzz[ix][iz]-tzz[ix-1][iz]) + C2*(tzz[ix+1][iz]-tzz[ix-2][iz]);
	  dtzx = C1*(txz[ix][iz]-txz[ix][iz-1]) + C2*(txz[ix][iz+1]-txz[ix][iz-2]);
	  vz[ix][iz] = vz[ix][iz] + ( dtx*dtzx + dtz*dtzz ) * rhoinv[ix][iz];
	 }
  }
  return;
}

void update_stress_iso(int nz, int nx, 
	                   float **vx, float **vz, float **txx, float **tzz, float **txz,
	                   float **c11, float **c55, float dtx, float dtz)
{
  int iz,ix;
  float dvxx,dvxz,dvzz,dvzx;
  for (iz=2;iz<nz-2;iz++) {
    for (ix=2;ix<nx-2;ix++) {
	  dvzz = C1*(vz[ix+1][iz]-vz[ix][iz]) + C2*(vz[ix+2][iz]-vz[ix-1][iz]);
	  dvxx = C1*(vx[ix][iz]-vx[ix][iz-1]) + C2*(vx[ix][iz+1]-vx[ix][iz-2]);
	  txx[ix][iz] = txx[ix][iz] + c11[ix][iz]*dtx*dvxx + (c11[ix][iz]-2*c55[ix][iz])*dtz*dvzz;
	  tzz[ix][iz] = tzz[ix][iz] + c11[ix][iz]*dtz*dvzz + (c11[ix][iz]-2*c55[ix][iz])*dtx*dvxx;
	  dvzx = C1*(vz[ix][iz+1]-vz[ix][iz]) + C2*(vz[ix][iz+2]-vz[ix][iz-1]);
	  dvxz = C1*(vx[ix][iz]-vx[ix-1][iz]) + C2*(vx[ix+1][iz]-vx[ix-2][iz]);
	  txz[ix][iz] = txz[ix][iz] + c55[ix][iz]*(dtz*dvxz + dtx*dvzx);
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

void v_to_p(float **p, float **vx, float **vy, float **vz, int nz, int nx, float dt, float dz, float dx, float vp_water, float rho_water)
{
  int iz, ix;

  for (iz=1;iz<nz-1;iz++){ 
    for (ix=1;ix<nx-1;ix++){
      /* here vy is left out since we are only considering the z-x plane (2d) */
	  p[ix][iz] = -dt*vp_water*vp_water*rho_water*(fd_approx_deriv1_order2(vz[ix][iz-1],vz[ix][iz+1],dz) + 
	  								               fd_approx_deriv1_order2(vx[ix-1][iz],vx[ix+1][iz],dx));
    }
  }
  for (iz=0;iz<nz;iz++){ 
	  p[0][iz]  = 0;
	  p[nx-1][iz] = 0;
  }
  for (ix=0;ix<nx;ix++){ 
	  p[ix][0]  = 0;
	  p[ix][nz-1] = 0;
  }
  return;
}

float fd_approx_deriv1_order2(float f1, float f3, float dx)
{
  /* 2nd order finite difference approximation to 1st derivative*/
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
