clear all;
tic;
c0=3e8;
delt = 1;
delx = 2*c0*delt;
dely = delx;
tsteps=2000;
xsteps=300;
ysteps=500;
eps_zz = ones(ysteps,xsteps);
mu_xx = ones(ysteps,xsteps);
mu_yy = mu_xx;
Ez=zeros(ysteps,xsteps);
Dz=zeros(ysteps,xsteps);
Hx=zeros(ysteps,xsteps);
Hy=zeros(ysteps,xsteps);
xsource = 50;
ysource = 50;


src_lambda = 100*delx;
src_omega = 2*pi*c0/src_lambda;
src_period = 0.2*src_lambda/c0;
src_freq = 1/src_period;
max_freq = 0.5*src_freq;


time = linspace(0,delt*tsteps, tsteps);
%Ey_source = sin(src_omega.*time);
nsrc = 1;
Ez_source = exp(-((time-6*src_period)/(src_period)).^2);

N=0;
Nspace = 30;

close all;
figure;

up_hx = -(c0*delt)./mu_xx;
up_hy = (c0*delt)./mu_yy;
up_dz = c0*delt*ones(ysteps,xsteps);
up_ez = 1./eps_zz;


for t=1:tsteps
	
%========update H field===========
	CEx(1:ysteps-1,1:xsteps)=diff(Ez,1,1)/dely;
% upper-y - boundary condition	
	CEx(ysteps,1:xsteps)= (0-Ez(ysteps,1:xsteps))/dely;

	CEy(1:ysteps, 1:xsteps-1) = diff(Ez,1,2)/delx;
% right-x boundary condition
	CEy(1:ysteps, xsteps) = (0 - Ez(1:ysteps,xsteps))/delx;
%update Hx	
	Hx+=up_hx.*CEx;
%update Hy
	Hy+=up_hy.*CEy;

%============update E field=========

	CHz(2:ysteps,2:xsteps)=diff(Hy(2:ysteps,:),1,2)/delx - diff(Hx(:,2:xsteps),1,1)/dely;
%boundary at low y
	CHz(1,2:xsteps)= diff(Hy(1,:))/delx - (Hx(1,2:xsteps)-0)/dely;
%boundaty at low x
	CHz(2:ysteps,1)= (Hy(2:ysteps,1)-0)/delx - diff(Hx(:,1))/dely;

%corner boundaty at (1, 1) 
	CHz(1, 1) = (Hy(1,1)-0)/delx - (Hx(1,1)-0)/dely;

%update Dz
	Dz+=up_dz.*CHz;
%update Ez
	Ez=up_ez.*Dz;

Ez(ysource, xsource)+=Ez_source(t);

%============plot the fields
	if N>Nspace
%		ylim([-1.5,1.5],"manual");		
		imshow(Ez,[-1e-3,1e-3]);
		title(sprintf('timestep is %d', t))
		pause(0.001);
		N = 0;		
	end;
	N++;
end
toc;


