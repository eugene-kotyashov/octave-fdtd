clear all;
tic;
c0 = physical_constant("speed of light in vacuum");
qe = physical_constant("elementary charge");
me = physical_constant("electron mass");
eps0 = physical_constant("electric constant");
hbar = physical_constant("Planck constant")/(2*pi);
eta0 = physical_constant("characteristic impedance of vacuum");
mu0 = physical_constant("mag. constant");

tsteps=200;
zsteps=500;
eps_yy = ones(1, zsteps);
mu_xx = ones(1, zsteps);
Ey=zeros(1,zsteps);
Dy=zeros(1,zsteps);
Py=zeros(1,zsteps);
delN = zeros(1,zsteps);
Py_old = zeros(1,zsteps);
delN_old = zeros(1,zsteps);

Hx=zeros(1,zsteps);
zsource = 50;
%device specification 
start_device = 100;
end_device = 300;
eps_device = 1;
eps_yy(start_device:end_device) = eps_device; 

%two level parameters
omega_a_val = 2*pi*5e14;
omega_a = zeros(1,zsteps);
omega_a(start_device:end_device) = omega_a_val;;
del_omega_a = zeros(1,zsteps);
del_omega_a(start_device:end_device) = 2*pi*5e13;
nk1 = zeros(1,zsteps);
nk1(start_device:end_device) = mu0./(hbar*omega_a(start_device:end_device));
delN0_val = 1e26;
delN0=zeros(1,zsteps);
delN0(start_device:end_device) = delN0_val;
tau21 = 20e-9;
invtau21 = zeros(1,zsteps);
invtau21(start_device:end_device) =1/tau21;
gamma_r=1e7;
gamma_CEO=(qe^2/me)*(omega_a.^2/(6*pi*eps0*c0^3));
kappa = zeros(1,zsteps);
kappa(start_device:end_device) = ((1/eps0)*gamma_r./gamma_CEO(start_device:end_device))*(qe^2/me);
%define simulation steps
delt = 0.005*(2*pi/omega_a_val);
delz = 2*c0*delt;

up_hx = ((c0*delt)/delz)./mu_xx;
up_dy = (c0*delt)/delz;

src_lambda = 1000*delz;
src_omega = 2*pi*c0/src_lambda;
src_period = 0.2*src_lambda/c0;
src_freq = 1/src_period;

max_freq = 2*omega_a_val/(2*pi);

%prepare fourier transform arrays
nfreqs = 400;
freqs = linspace(0,max_freq,nfreqs);
K = exp(-1i*2*pi*delt*freqs);
ref = zeros(1, nfreqs);
trans = zeros(1, nfreqs);
norm_src = zeros(1, nfreqs);


time = linspace(0,delt*tsteps, tsteps);
nsrc = 1;
Ey0=1e7;
Ey_source = Ey0*exp(-((time-6*src_period)/(src_period)).^2).*sin(omega_a_val*time);
Hx_source = -Ey0*sqrt(eps_yy(1,zsource)/mu_xx(1,zsource)).*exp(-(((time+0.5*delz*nsrc/c0+0.5*delt)-6*src_period)/(src_period)).^2).*sin(omega_a_val*time);
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
%two level additional equations for polatization and inversion
	Py_old = Py;	
	Py = ((del_omega_a*delt-2).*Py_old + 2*delt^2*kappa.*delN.*Ey + (4-2*(omega_a*delt).^2).*Py)./(del_omega_a*delt+2);
	delN_old = delN;
	delN= - 2*nk1.*Ey.*(Py-Py_old) + delN_old + 2*delt*(delN0-delN).*invtau21;
			

%       ey boundary condition
	Dy(1,1)+=up_dy*(Hx(1,1)-h3);

	Dy(2:zsteps)+=up_dy*diff(Hx);

	e3 = e2;
	e2 = e1;
	e1 = Dy(1,zsteps);
%	TF/SF source correction
	Dy(1,zsource)-=up_dy*Hx_source(1,t);
%	soft source
%	Ey(1,zsource)+=Ey_source(1,t);

%	calculation Ey from Dy
	%Ey = Dy./eps_yy;
	Ey = (Dy-Py)./eps_yy;

	tmp = K.^t;
	ref+= Ey(1)*tmp;
	trans+=Ey(zsteps)*tmp;
	norm_src+=Ey_source(t)*tmp;		


	if N>Nspace
		subplot(3,1,1);		
		plot( linspace(1,zsteps,zsteps),Ey,"r",
			linspace(1,zsteps,zsteps),Hx,"b");
		title(sprintf("%d fields", t));
		legend("Ey", "Hz");		
		ylim([-1.5*Ey0,1.5*Ey0],"manual");
		subplot(3,1,2);
		plot( linspace(1,zsteps,zsteps),delN,"r");%-delN0)/delN0_val,"r");
		title("population difference");		
		subplot(3,1,3);
		plot(2*pi*freqs,(abs(ref)./abs(norm_src)).^2,"m",
			2*pi*freqs,(abs(trans)./abs(norm_src)).^2,"g",
			2*pi*freqs, (abs(ref)./abs(norm_src)).^2 + (abs(trans)./abs(norm_src)).^2,"b");		
		title("reflectance and transmittance");
		legend("Ref", "Trans", "Tot");
		ylim([0,1.2],"manual");
		pause(0.001);
		N = 0;		
	end;
	N++;
end
toc;
%ref=ref*delt;
%trans=trans*delt;
%norm_src=norm_src*delt;


