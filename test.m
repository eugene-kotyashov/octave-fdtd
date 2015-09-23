clear all;
tic;
c0=3e8;
delt = 1;
delz = 2*c0*delt;
tsteps=10000;
zsteps=500;
eps_yy = ones(1, zsteps);
mu_xx = ones(1, zsteps);
Ey=zeros(1,zsteps);
Hx=zeros(1,zsteps);
zsource = 50;
%device specification 
start_device = 100;
end_device = 300;
eps_device = 9;
eps_yy(start_device:end_device) = eps_device; 

up_hx = ((c0*delt)/delz)./eps_yy;
up_ey = ((c0*delt)/delz)./mu_xx;

src_lambda = 100*delz;
src_omega = 2*pi*c0/src_lambda;
src_period = 0.2*src_lambda/c0;
src_freq = 1/src_period;
max_freq = 0.5*src_freq;

%prepare fourier transform arrays
nfreqs = 400;
freqs = linspace(0,max_freq,nfreqs);
K = exp(-1i*2*pi*delt*freqs);
ref = zeros(1, nfreqs);
trans = zeros(1, nfreqs);
norm_src = zeros(1, nfreqs);


time = linspace(0,delt*tsteps, tsteps);
nsrc = 1;
Ey_source = exp(-((time-6*src_period)/(src_period)).^2);
Hx_source = -sqrt(eps_yy(1,zsource)/mu_xx(1,zsource)).*exp(-(((time+0.5*delz*nsrc/c0+0.5*delt)-6*src_period)/(src_period)).^2);
N=0;
Nspace = 100;

close all;
figure;
hax = gca();

h1 = h2 = h3 = 0;
e1 = e2 = e3 = 0;

for t=1:tsteps-1
	Hx(1:zsteps-1)+=up_hx(1:zsteps-1).*diff(Ey);			

%	hx boundary condition
	Hx(1,zsteps)+=up_hx(1,zsteps).*(e3-Ey(1,zsteps));
	h3 = h2;
	h2 = h1;
	h1 = Hx(1,1);
%	TF/SF source correction
	Hx(1,zsource-1)-=up_hx(1,zsource-1)*Ey_source(1,t);
%       ey boundary condition
	Ey(1,1)+=up_ey(1,1).*(Hx(1,1)-h3);

	Ey(2:zsteps)+=up_ey(2:zsteps).*diff(Hx);

	e3 = e2;
	e2 = e1;
	e1 = Ey(1,zsteps);
%	TF/SF source correction
	Ey(1,zsource)-=up_ey(1,zsource)*Hx_source(1,t);
%	soft source
%	Ey(1,zsource)+=Ey_source(1,t);

	tmp = K.^t;
	ref+= Ey(1)*tmp;
	trans+=Ey(zsteps)*tmp;
	norm_src+=Ey_source(t)*tmp;		


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
end
toc;
%ref=ref*delt;
%trans=trans*delt;
%norm_src=norm_src*delt;


