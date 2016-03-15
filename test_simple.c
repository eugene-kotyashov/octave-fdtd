#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <pthread.h>







pthread_barrier_t time_step_barrier;
int NUM_THREADS = 2;

typedef struct tag_chunk_info_t
{
	int 	tid;	
} chunk_info_t;



double * Ey;
double * up_hx;
double * up_ey;
double * Hx;
double * Ey_source;
double * Hx_source;
int	tsteps;
int	zsteps;
int	zsource;

double *Ey_source;
double *Ey_time;

double *Hx_source;

int nfreqs;
complex double * K;
complex double * ref;
complex double * trans;
complex double * norm_src;

double N=0;
double Nspace = 100;


double h1 = 0;
double h2 = 0;
double h3 = 0;
double e1 = 0;
double e2 = 0;
double e3 = 0;


#ifdef GNUPLOT_PIPING
FILE * gnuplotPipe;
#endif


complex double c_pow(complex double val, int power)
{
	double resRe, resIm;
	double rho =  pow(cabs(val),power);
	double phi = atan2(cimag(val),creal(val));//carg(val);
	resRe = rho*cos(phi*power);
	resIm = rho*sin(phi*power);
	return resRe + I*resIm;
}

complex double c_exp(complex double val)
{
	double resRe, resIm;
	double valReExp = exp(creal(val));
	double valIm = cimag(val);
	resRe = valReExp*cos(valIm);
	resIm = valReExp*sin(valIm);
	return resRe + I*resIm;

}


void * do_time_step(void * input)
{

chunk_info_t * info = (chunk_info_t*)input;
if (info == NULL) pthread_exit(NULL);

int freq_start = nfreqs*info->tid/NUM_THREADS;
int freq_end = nfreqs*(info->tid+1)/NUM_THREADS-1;

//for Hx update iz runs from 0 to zsteps-2 - total of zsteps -1 elements to update
int z_start_hx = (zsteps - 1)*info->tid/NUM_THREADS;
int z_end_hx = (zsteps - 1)*(info->tid+1)/NUM_THREADS-1;	
//for Ey update iz runs from 1 to zsteps-1 - total of zsteps -1 elements to update
int z_start_ey = z_start_hx + 1;
int z_end_ey = z_end_hx + 1;

int do_ey_source = (zsource <= z_end_ey)&&(zsource >= z_start_ey);
int do_hx_source = ( (zsource -1) <= z_end_hx)&&( (zsource-1) >= z_start_hx);

if (info->tid == NUM_THREADS-1)
{
	z_end_hx = zsteps-2;
	z_end_ey = zsteps-1;
}

int t;
for ( t=0; t<tsteps; t++)
{
	
	int iz;
//========update Hx========================
	for(iz=z_start_hx;iz<=z_end_hx;iz++)
	{	
		Hx[iz]+=up_hx[iz]*(Ey[iz+1]-Ey[iz]);			
	}
//	hx boundary condition - right boundary
	if (z_end_hx == zsteps-2)
	{
		Hx[zsteps-1]+=up_hx[zsteps-1]*(e3-Ey[zsteps-1]);
	}
//collect data needed for ey boundary condition
	if (z_start_hx == 0)
	{
		h3 = h2;
		h2 = h1;
		h1 = Hx[0];
	}
//	TF/SF source correction		
	if (do_hx_source)
	{
		Hx[zsource-1] -= up_hx[zsource-1]*Ey_source[t];
	}
//wait for all threads done with Hx
	pthread_barrier_wait(&time_step_barrier);

//===========update Ey=====================
//      ey boundary condition
	if (z_start_ey == 1)
	{
		Ey[0]+=up_ey[0]*(Hx[0]-h3);
	}

	for(iz=z_start_ey;iz<=z_end_ey;iz++)
	{
		Ey[iz]+=up_ey[iz]*(Hx[iz]-Hx[iz-1]);
	}
//collect data needed for ey boundary condition
	if(z_end_ey == zsteps-1)
	{
		e3 = e2;
		e2 = e1;
		e1 = Ey[zsteps-1];
	}
//	TF/SF source correction
	if (do_ey_source)
	{
		Ey[zsource]-=up_ey[zsource]*Hx_source[t];
	}
//wait all threads done wuth Ey
	pthread_barrier_wait(&time_step_barrier);

//perform Fourier transform

	complex double tmp;
	for(iz=freq_start;iz<=freq_end; iz++){
		tmp = c_pow(K[iz],t);
		ref[iz]+= Ey[0]*tmp;
		trans[iz]+=Ey[zsteps-1]*tmp;
		norm_src[iz]+=Ey_source[t]*tmp;		
	}
	if (info->tid == NUM_THREADS-1)
		Ey_time[t] = Ey[zsteps-1];

#ifdef GNUPLOT_PIPING
	if ( (N>Nspace)&&(info->tid == 0) )
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

pthread_exit(NULL);

}


int main(int argc, char* args[])
{

if (argc == 2) 
{
	sscanf(args[1],"%d",&NUM_THREADS);
}

printf("\nusing %d threads\n", NUM_THREADS);

double pi = M_PI;
double c0=3e8;
double delt = 1;
double delz = 2*c0*delt;

tsteps=200000;
zsteps=500;

double* eps_yy = (double*)malloc(zsteps*sizeof(double));
double* mu_xx = (double*)malloc(zsteps*sizeof(double));
Ey=(double*)malloc(zsteps*sizeof(double));
Hx=(double*)malloc(zsteps*sizeof(double));
up_hx = (double*)malloc(zsteps*sizeof(double));
up_ey = (double*)malloc(zsteps*sizeof(double));
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
zsource = 50;
double src_lambda = 100*delz;
double src_omega = 2*pi*c0/src_lambda;
double src_period = 0.2*src_lambda/c0;
double src_freq = 1/src_period;
double max_freq = 0.5*src_freq;

//prepare fourier transform arrays
nfreqs = 400;
K = (complex double*)malloc(nfreqs*sizeof(complex double));
ref = (complex double*)malloc(nfreqs*sizeof(complex double));
trans = (complex double*)malloc(nfreqs*sizeof(complex double));
norm_src = (complex double*)malloc(nfreqs*sizeof(complex double));

double del_freq = max_freq/nfreqs;
for(iz=0;iz<nfreqs; iz++)
{
	K[iz] = c_exp(-I*2*pi*delt*del_freq*iz);
	ref[iz] = 0;
	trans[iz] = 0;
	norm_src[iz] = 0;
}


Ey_source = (double*)malloc(tsteps*sizeof(double));
Ey_time = (double*)malloc(tsteps*sizeof(double));

Hx_source = (double*)malloc(tsteps*sizeof(double));
double time = 0;
double nsrc = sqrt(eps_yy[zsource]*mu_xx[zsource]);

for(iz=0;iz<tsteps;iz++)
{
	time = delt*iz;	
	Ey_source[iz] = exp(-((time-6*src_period)/(src_period))*((time-6*src_period)/(src_period)));
	Hx_source[iz] = -sqrt(eps_yy[zsource]/mu_xx[zsource])*exp(-((time+0.5*delz*nsrc/c0+0.5*delt)-6*src_period)/(src_period)*((time+0.5*delz*nsrc/c0+0.5*delt)-6*src_period)/(src_period));
}



#ifdef GNUPLOT_PIPING
gnuplotPipe = popen ("gnuplot ", "w");
#endif

int t;

pthread_t threads[NUM_THREADS];
pthread_attr_t attr;
pthread_attr_init(&attr);
pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
int rc;
void* status;
chunk_info_t info[NUM_THREADS];

pthread_barrier_init(&time_step_barrier,NULL,NUM_THREADS);

for(iz=0;iz<NUM_THREADS;iz++) 
{
	info[iz].tid = iz;
	rc = pthread_create(&threads[iz], &attr, do_time_step, (void*)&info[iz]);
	if (rc) 
	{
	        printf("ERROR; return code from pthread_create() is %d\n", rc);
	        exit(-1);
        }
}


for(iz=0;iz<NUM_THREADS;iz++) 
{
	info[iz].tid = iz;
	rc = pthread_join(threads[iz], &status);
	if (rc) {
	        printf("ERROR; return code from pthread_join() is %d\n", rc);
        exit(-1);
        }
}



pthread_attr_destroy(&attr);
pthread_barrier_destroy(&time_step_barrier);




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
