clear all;
tic;
c0 = physical_constant("speed of light in vacuum");
qe = physical_constant("elementary charge");
me = physical_constant("electron mass");
eps0 = physical_constant("electric constant");
hbar = physical_constant("Planck constant")/(2*pi);
eta0 = physical_constant("characteristic impedance of vacuum");
mu0 = physical_constant("mag. constant");


%two level parameters
omega_a_val = 2*pi*5e14;



%define simulation steps
delt = 0.05*(2*pi/(3*omega_a_val));
eps_device = 4;
delz = 2*delt*c0;

tsteps=1000000;
%device specification 
start_device = 100;
Length = 0.3e-5;
end_device = start_device+floor(Length/delz);
zsteps=end_device + 50;

eps_yy = ones(1, zsteps);
mu_xx = ones(1, zsteps);
Ey=zeros(1,zsteps);
Dy=zeros(1,zsteps);
Py=zeros(1,zsteps);
delN = zeros(1,zsteps);
Py_old2 = zeros(1,zsteps);
Py_old = zeros(1,zsteps);
delN_old = zeros(1,zsteps);
delN_old2 = zeros(1,zsteps);
Eyold = zeros(1,zsteps);
Hx=zeros(1,zsteps);
zsource = 50;


eps_yy(start_device:end_device) = eps_device; 

%two level parameters
omega_a = zeros(1,zsteps);
omega_a(start_device:end_device) = omega_a_val;
del_omega_a_val = 2*pi*5e13;
del_omega_a = zeros(1,zsteps);
del_omega_a(start_device:end_device) = del_omega_a_val;
nk1 = zeros(1,zsteps);
nk1(start_device:end_device) = mu0./(hbar*omega_a(start_device:end_device));
delN0_val = 1e26;
delN0=zeros(1,zsteps);
delN0(start_device:end_device) = delN0_val;
delN = -delN0;
gamma_r=5e7; %radiative decay rate
tau21 = 1/gamma_r;
invtau21 = zeros(1,zsteps);
invtau21(start_device:end_device) =1/tau21;

gamma_CEO=(qe^2/me)*(omega_a.^2/(6*pi*eps0*c0^3));
gamma_CEO_val=(qe^2/me)*(omega_a_val^2/(6*pi*eps0*c0^3));
kappa = zeros(1,zsteps);
kappa_val = ((1/eps0)*gamma_r/gamma_CEO_val)*(qe^2/me);
kappa(start_device:end_device) = ((1/eps0)*gamma_r./gamma_CEO(start_device:end_device))*(qe^2/me);
plasma_freq_val = 2*pi*sqrt(kappa_val*delN0_val);

%update coeffs for delN
n1 = zeros(1,zsteps);
n1_val = (delt-2*tau21)/(delt+2*tau21);
n1(start_device:end_device)=n1_val;

n2 = zeros(1,zsteps);
n2_val = 2*tau21*eta0/((c0*omega_a_val*hbar)*(delt+2*tau21));
n2(start_device:end_device)=n2_val;

n3 = zeros(1,zsteps);
n3_val = 2*delN0_val*delt/(delt+2*tau21);
n3(start_device:end_device)=n3_val;

pump_rate = 10000*delN0_val/tau21;
n4 = zeros(1,zsteps);
n4_val = 2*pump_rate*tau21*delt/(delt+2*tau21);
n4(start_device:end_device)=n4_val;


up_hx = ((c0*delt)/delz)./mu_xx;
up_dy = (c0*delt)/delz;



max_freq = 3*omega_a_val/(2*pi);

%prepare fourier transform arrays
nfreqs = 600;
freqs = linspace(0,max_freq,nfreqs);
omegas = 2*pi*freqs;
K = exp(-1i*2*pi*delt*freqs);
ref = zeros(1, nfreqs);
trans = zeros(1, nfreqs);
norm_src = zeros(1, nfreqs);
delN_ft = zeros(1, nfreqs);
delN_time =zeros(1,tsteps);
_time =zeros(1,tsteps);


time_t = linspace(0,delt*tsteps, tsteps);
nsrc = 1;
src_omega = omega_a_val;
src_period = 2*pi/src_omega;

Ey0=1e3;
Ey_source = Ey0*exp(-((time_t-6*src_period)/(src_period)).^2).*sin(omega_a_val*time_t);
Hx_source = -Ey0*sqrt(eps_yy(1,zsource)/mu_xx(1,zsource)).*exp(-(((time_t+0.5*delz*nsrc/c0+0.5*delt)-6*src_period)/(src_period)).^2).*sin(omega_a_val*(time_t+0.5*delz*nsrc/c0+0.5*delt));
N=0;
Nspace = 100;
delN_time_t = zeros(1,tsteps);

close all;
%figure;

plasma_freq_val = sqrt(kappa_val*delN0_val);


% ------------------- calculation of reflection coefficients
%eps_w = zeros(1,nfreqs)+eps_device;
eps_w = eps_device+plasma_freq_val^2./(omega_a_val^2 - omegas.^2 - I*del_omega_a_val*omegas);
ref_w = zeros(1,nfreqs);
trans_w = zeros(1,nfreqs);
e_src=[0; 1];
ey_trans = [0; 0];
ey_ref = [0; 0];
k_norm = omega_a_val/c0;
%L = 2*(2*pi/k_norm);
L = (end_device - start_device)*delz;
S0_11 = zeros(2);
S0_12 = eye(2);
S0_21 = S0_12;
S0_22 = S0_11;
for i_w = 1:nfreqs
	k_0 = omegas(i_w)/c0;
	k_z = sqrt(eps_w(i_w));
	Q=[0, eps_w(i_w); -eps_w(i_w),0];
	OMEGA = I*k_z*eye(2);
	Qh = [0 1; -1 0];
	Vh = -I*Qh;
	V = Q*(OMEGA^(-1));
	X=expm(OMEGA*k_0*L);
	A=eye(2)+V^(-1)*Vh;
	B=eye(2)-V^(-1)*Vh;
	D=A-X*B*(A^(-1))*X*B;
	S11=(D^(-1))*(X*B*(A^(-1))*X*A - B);
	S22=S11;
	S12=(D^(-1))*X*(A-B*(A^(-1))*B);
	S21=S12;
%	S = redheffer(S0_11,S0_12,S0_21,S0_22,S11,S12,S21,S22);
%	S11=[S(1,1) S(1,2); S(2,1) S(2,2)];
%	S12=[S(3,1) S(3,2); S(4,1) S(4,2)];
%	S21=[S(5,1) S(5,2); S(6,1) S(6,2)];
%	S22=[S(7,1), S(7,2); S(8,1) S(8,2)];
	ey_ref = S11*e_src;
	ey_trans = S21*e_src;
	ref_w(i_w)=abs(ey_ref(2))^2;
	trans_w(i_w)=abs(ey_trans(2))^2;
	
end
%-------------------------------------------


h1 = h2 = h3 = 0;
e1 = e2 = e3 = 0;


for t=1:tsteps
	Hx(1:zsteps-1)+=up_hx(1:zsteps-1).*diff(Ey);			

%	hx boundary condition
	Hx(1,zsteps)+=up_hx(1,zsteps).*(e3-Ey(1,zsteps));
	h3 = h2;
	h2 = h1;
	h1 = Hx(1,1);
%	TF/SF source correction
	Hx(1,zsource-1)-=up_hx(1,zsource-1)*Ey_source(1,t);
%two level additional equations for polatization and inversion
	Py_old2 = Py_old;	
	Py_old = Py;
	Py = ((del_omega_a*delt-2).*Py_old2 + 2*delt^2*kappa.*delN.*Ey + (4-2*(omega_a*delt).^2).*Py)./(del_omega_a*delt+2);
	%delN_old2 = delN_old;
	%delN_old = delN;
	%delN= - 2*nk1.*Ey.*(Py-Py_old2) + delN_old2 + 2*delt*(delN0-delN).*invtau21;
			

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
	Eyold = Ey;
	Ey = (Dy-Py)./eps_yy;

%	update delN here
	delN = - n1.*delN - n2.*(Ey+Eyold).*(Py-Py_old) + n3 - n4;

	tmp = K.^t;
	ref+= Ey(1)*tmp;
	trans+=Ey(zsteps)*tmp;
	norm_src+=Ey_source(t)*tmp;
	if t>1000
		delN_ft+=delN(ceil(0.5*(start_device+end_device)))*tmp/delN0_val;
	end;
	delN_time(t) = delN(ceil(0.5*(start_device+end_device)));
	Ey_time(t)=Ey(zsteps);


%	if N==Nspace
%		subplot(3,1,1);		
%		plot( linspace(1,zsteps,zsteps),Ey,"r",
%			linspace(1,zsteps,zsteps),Hx,"b");
%		title(sprintf("%d fields", t));
%		legend("Ey", "Hz");		
%		ylim([-1.5*Ey0,1.5*Ey0],"manual");
%		subplot(3,1,2);
%		plot( linspace(1,zsteps,zsteps),delN,"r");
%		title("population difference in time");		
%		subplot(3,1,3);
%		plot(
%			2*pi*freqs,(abs(ref)./abs(norm_src)).^2,"m",
%			2*pi*freqs,(abs(trans)./abs(norm_src)).^2,"g",
%			2*pi*freqs, (abs(ref)./abs(norm_src)).^2 + (abs(trans)./abs(norm_src)).^2,"b",
%			omegas,trans_w,"c-"
%		);		
%		title("reflectance and transmittance");
%		legend("Ref", "Trans", "Tot", "calc_trans");
%		ylim([0,1.2],"manual");
%		pause(0.001);		
%		N = 0;
%	end;
%	N++;
		delN_time(t) = delN(ceil(0.5*(start_device+end_device)));	
end
plot(
	2*pi*freqs,(abs(ref)./abs(norm_src)).^2,"m",
	2*pi*freqs,(abs(trans)./abs(norm_src)).^2,"g",
	2*pi*freqs, (abs(ref)./abs(norm_src)).^2 + (abs(trans)./abs(norm_src)).^2,"b",
	omegas,trans_w,"c-"
);		
title("reflectance and transmittance");
legend("Ref", "Trans", "Tot", "Ctrans");
ylim([0,1.2],"manual");
figure;
subplot(2,1,1);
plot( time_t,delN_time,"r");
title("population difference in time");
subplot(2,1,2);
plot( time_t,Ey_time,"b");
title("Ey transmitterd in time");		
toc;
%ref=ref*delt;
%trans=trans*delt;
%norm_src=norm_src*delt;





