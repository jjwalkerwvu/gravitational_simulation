% particle_init.m

% this function initializes the particles for the simulation.


% try first: homogeneous spatial distribution, maxwellian velocity distribution

% for spatial distribution, try random distribution in some spherical volume?

function [x0,y0,z0,vx0,vy0,vz0]=particle_init(x,y,z,vx,vz,vy,lsim)

% find the number of particles in the simulation volume:
nparticles=length(x);



x0=lsim*rand(size(x));
y0=lsim*rand(size(x));
z0=lsim*rand(size(x));
vx0=zeros(size(x));
vy0=zeros(size(x));
vz0=zeros(size(x));

% thermal speed? sqrt(2*qe*T/m)?
vmax=1;
% quick something I tried that sort of works:
xplot=linspace(-vmax,vmax,nparticles);
% to eliminate a lot of the tail, just reduce by some factor?
vx0=erf(xplot);
%figure(1);clf;plot(vx0,exp(-vx0.^2));
%vx0=zeros(size(x));

%f=linspace(0,1,length(x));
%vx0=1/2*(1+erf(xplot/sqrt(2)));
vy0=vx0;
vz0=vx0;