#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

void main()
{
double pi = M_PI;
double c0=3e8;
double delt = 1;
double delz = 2*c0*delt;
int tsteps=100000;
int zsteps=500;

double* eps_yy = (double*)malloc(zsteps*sizeof(double));
double* mu_xx = (double*)malloc(zsteps*sizeof(double));
double* Ey=(double*)malloc(zsteps*sizeof(double));
double* Hx=(double*)malloc(zsteps*sizeof(double));
double* up_hx = (double*)malloc(zsteps*sizeof(double));
double* up_ey = (double*)malloc(zsteps*sizeof(double));
double k_grid = (c0*delt)/delz;




//device specification 
int start_device = 100;
int end_device = 300;
int eps_device = 9;
int iz = 0;

for (iz=0;iz<zsteps;iz++)
{
	Ey[iz] = 0;
	Hx[iz] = 0;
	eps_yy[iz] = 1.0;
	mu_xx[iz] = 1.0;
	up_hx[iz] = k_grid;
	up_ey[iz] = k_grid;	
}
//set device parameters
for (iz=start_device;iz<=end_device;iz++)
{
	eps_yy[iz] = eps_device;	
}
//set update coefficients
for(iz=0;iz<zsteps;iz++)
{
	up_hx[iz] /= eps_yy[iz];
	up_ey[iz] /= mu_xx[iz];
}

//calculate source
int zsource = 50;
double src_lambda = 100*delz;
double src_omega = 2*pi*c0/src_lambda;
double src_period = 0.2*src_lambda/c0;
double src_freq = 1/src_period;
double max_freq = 0.5*src_freq;

//prepare fourier transform arrays
int nfreqs = 400;
double complex * K = (double complex*)malloc(nfreqs*sizeof(double complex));
double complex * ref = (double complex*)malloc(nfreqs*sizeof(double complex));
double complex * trans = (double complex*)malloc(nfreqs*sizeof(double complex));
double complex * norm_src = (double complex*)malloc(nfreqs*sizeof(double complex));
double del_freq = max_freq/nfreqs;
for(iz=0;iz<nfreqs; iz++)
{
	K[iz] = exp(I*2*pi*delt*del_freq*iz);
	ref[iz] = 0;
	trans[iz] = 0;
	norm_src[iz] = 0;
}


double *Ey_source = (double*)malloc(tsteps*sizeof(double));

double *Hx_source = (double*)malloc(tsteps*sizeof(double));
double time = 0;
double nsrc = sqrt(eps_yy[zsource]*mu_xx[zsource]);

for(iz=0;iz<tsteps;iz++)
{
	time = delt*iz;	
	Ey_source[iz] = exp(-((time-6*src_period)/(src_period))*((time-6*src_period)/(src_period)));
	Hx_source[iz] = -sqrt(eps_yy[zsource]/mu_xx[zsource])*exp(-((time+0.5*delz*nsrc/c0+0.5*delt)-6*src_period)/(src_period)*((time+0.5*delz*nsrc/c0+0.5*delt)-6*src_period)/(src_period));
}

double N=0;
double Nspace = 100;


double h1 = 0;
double h2 = 0;
double h3 = 0;
double e1 = 0;
double e2 = 0;
double e3 = 0;

int t;
for ( t=0; t<tsteps; t++)
{

	
	for(iz=0;iz<zsteps-1;iz++)
	{	
		Hx[iz]+=up_hx[iz]*Ey[iz+1]-Ey[iz];			
	}

//	hx boundary condition
	Hx[zsteps-1]+=up_hx[zsteps-1]*(e3-Ey[zsteps-1]);
	h3 = h2;
	h2 = h1;
	h1 = Hx[0];
//	TF/SF source correction
	Hx[zsource-1]-=up_hx[zsource-1]*Ey_source[0];
//      ey boundary condition
	Ey[0]+=up_ey[0]*(Hx[0]-h3);

	for(iz=1;iz<zsteps;iz++)
	{
		Ey[iz]+=up_ey[iz]*(Hx[iz]-Hx[iz-1]);
	}

	e3 = e2;
	e2 = e1;
	e1 = Ey[0];
//	TF/SF source correction
	Ey[zsource]-=up_ey[zsource]*Hx_source[t];
//	soft source
//	Ey(1,zsource)+=Ey_source(1,t);
	double tmp = 0;
	for(iz=0;iz<nfreqs; iz++){
		tmp = pow(K[iz],t);
		ref[iz]+= Ey[0]*tmp;
		trans[iz]+=Ey[zsteps-1]*tmp;
		norm_src[iz]+=Ey_source[t]*tmp;		
	}
/*

	if N>Nspace
		subplot(2,1,1);		
		plot( linspace(1,zsteps,zsteps),Ey,"r",
			linspace(1,zsteps,zsteps),Hx,"b");
		title(num2str(t));		
		ylim([-1.5,1.5],"manual");
		subplot(2,1,2);
		plot(freqs,(abs(ref)./abs(norm_src)).^2,"m",
			freqs,(abs(trans)./abs(norm_src)).^2,"g",
			freqs, (abs(ref)./abs(norm_src)).^2 + (abs(trans)./abs(norm_src)).^2,"b");		
		title("reflectance and transmittance");		
		pause(0.001);
		N = 0;		
	end;
	N++;
}


*/
}
free(eps_yy);
free(mu_xx);
free(Ey);
free(Hx);
free(up_hx);
free(up_ey);

free(K);
free(ref);
free(trans);
free(norm_src);

free(Hx_source);
free(Ey_source);

}
