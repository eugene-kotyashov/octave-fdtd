#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <pthread.h>


#define NUM_THREADS 4


typedef struct tag_chunk_info_t
{
	int 	tid;	
	double * Ey;
	double * up_hx;
	double * up_ey;
	double * Hx;
	int 	z_start;
	int	z_end;
	
} chunk_info_t;

void* updateHx(void* input)
{		
	chunk_info_t * info = (chunk_info_t*)input;
	if (info == NULL) pthread_exit(NULL);
	int iz;
	for(iz=info->z_start;iz<=info->z_end;iz++)
	{	
		info->Hx[iz]+=info->up_hx[iz]*(info->Ey[iz+1]-info->Ey[iz]);			
	}
	pthread_exit(NULL);
}

void* updateEy(void* input)
{
	chunk_info_t * info = (chunk_info_t*)input;
	if (info == NULL) pthread_exit(NULL);
	int iz;
	for(iz=info->z_start;iz<=info->z_end;iz++)
	{
		info->Ey[iz]+=info->up_ey[iz]*(info->Hx[iz]-info->Hx[iz-1]);
	}
	pthread_exit(NULL);
}


void main()
{


double pi = M_PI;
double c0=3e8;
double delt = 1;
double delz = 2*c0*delt;
int tsteps=20000;
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
complex double * K = (complex double*)malloc(nfreqs*sizeof(complex double));
complex double * ref = (complex double*)malloc(nfreqs*sizeof(complex double));
complex double * trans = (complex double*)malloc(nfreqs*sizeof(complex double));
complex double * norm_src = (complex double*)malloc(nfreqs*sizeof(complex double));
double del_freq = max_freq/nfreqs;
for(iz=0;iz<nfreqs; iz++)
{
	K[iz] = cexp(-I*2*pi*delt*del_freq*iz);
	ref[iz] = 0;
	trans[iz] = 0;
	norm_src[iz] = 0;
}


double *Ey_source = (double*)malloc(tsteps*sizeof(double));
double *Ey_time = (double*)malloc(tsteps*sizeof(double));

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

#ifdef GNUPLOT_PIPING
FILE * gnuplotPipe = popen ("gnuplot ", "w");
#endif

int t;

pthread_t threads[NUM_THREADS];
pthread_attr_t attr;
pthread_attr_init(&attr);
pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
int rc;
void* status;
chunk_info_t hx_info[NUM_THREADS];
chunk_info_t ey_info[NUM_THREADS];
for(iz=0;iz<NUM_THREADS;iz++) 
{
	hx_info[iz].Ey = Ey;
	ey_info[iz].Ey = Ey;
	hx_info[iz].Hx = Hx;
	ey_info[iz].Hx = Hx;
	hx_info[iz].up_ey = up_ey;
	ey_info[iz].up_ey = up_ey;
	hx_info[iz].up_hx = up_hx;
	ey_info[iz].up_hx = up_hx;
	hx_info[iz].tid = iz;
//for Hx update iz runs from 0 to zsteps-2 - total of zsteps -1 elements to update
	hx_info[iz].z_start = (zsteps - 1)*iz/NUM_THREADS;
	hx_info[iz].z_end = (zsteps - 1)*(iz+1)/NUM_THREADS-1;	
//for Ey update iz runs from 1 to zsteps-1 - total of zsteps -1 elements to update
	ey_info[iz].z_start = (zsteps - 1)*iz/NUM_THREADS+1;
	ey_info[iz].z_end = (zsteps - 1)*(iz+1)/NUM_THREADS-1+1;
	if (iz == NUM_THREADS-1)
	{
		hx_info[iz].z_end = zsteps-2;
		ey_info[iz].z_end = zsteps-1;
	}
}

for ( t=0; t<tsteps; t++)
{

	for(iz=0;iz<NUM_THREADS;iz++)
	{		
		rc = pthread_create(&threads[iz],&attr, updateHx, (void*)(&hx_info[iz]));
	}
	for(iz=0;iz<NUM_THREADS;iz++)
	{
		rc = pthread_join(threads[iz],&status);
	}
	

//	hx boundary condition
	Hx[zsteps-1]+=up_hx[zsteps-1]*(e3-Ey[zsteps-1]);
	h3 = h2;
	h2 = h1;
	h1 = Hx[0];
//	TF/SF source correction
	Hx[zsource-1]-=up_hx[zsource-1]*Ey_source[t];
//      ey boundary condition
	Ey[0]+=up_ey[0]*(Hx[0]-h3);

	for(iz=0;iz<NUM_THREADS;iz++)
	{		
		rc = pthread_create(&threads[iz],&attr, updateEy, (void*)(&ey_info[iz]));
	}	
	for(iz=0;iz<NUM_THREADS;iz++)
	{
		rc = pthread_join(threads[iz],&status);
	}
	

	e3 = e2;
	e2 = e1;
	e1 = Ey[zsteps-1];
//	TF/SF source correction
	Ey[zsource]-=up_ey[zsource]*Hx_source[t];
//	soft source
//	Ey(1,zsource)+=Ey_source(1,t);
	complex double tmp;
	for(iz=0;iz<nfreqs; iz++){
		tmp = cpow(K[iz],t);
		ref[iz]+= Ey[0]*tmp;
		trans[iz]+=Ey[zsteps-1]*tmp;
		norm_src[iz]+=Ey_source[t]*tmp;		
	}
	Ey_time[t] = Ey[zsteps-1];

#ifdef GNUPLOT_PIPING
	if (N>Nspace)
	{			
		fprintf(gnuplotPipe,"plot '-' with lines\n");
		for(iz=0;iz<zsteps;iz++)
		{
			fprintf(gnuplotPipe, "%d %e\n",
				iz,
				Ey[iz]
			);
		}
		fprintf(gnuplotPipe,"e\n");
		N = 0;		
		usleep(1000000);
		
	}
	N++;
#endif

}

FILE* fout = fopen("out.txt","w");
for(iz=0;iz<nfreqs; iz++)
{
	fprintf(fout,"%e %e %e %e\n",
		2*pi*del_freq*iz,
		cabs(norm_src[iz]),
		cabs(ref[iz])*cabs(ref[iz])/(cabs(norm_src[iz])*cabs(norm_src[iz])),
		cabs(trans[iz])*cabs(trans[iz])/(cabs(norm_src[iz])*cabs(norm_src[iz]))
	);
}
fclose(fout);

fout = fopen("time_out.txt", "w");
for(iz=0;iz<tsteps;iz++)
{
	fprintf(fout,"%e %e %e\n",
		iz*delt,
		Ey_time[iz],
		Ey_source[iz]
	);	
}
fclose(fout);

#ifdef GNUPLOT_PIPING
pclose(gnuplotPipe);
#endif

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

pthread_exit(NULL);

}
