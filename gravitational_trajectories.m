% gravitational_trajectories.m

% this is the main program.
% using some input file? this script produces particle trajectories
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% maybe make x,y,z, etc. values 2d arrays? 

nparticles=1e1;
% bulk density of the particle material, set to one for now but might want to
% think more on what to set this to. Maybe set equal to water density in 
% whatever units?
rho=1;
nsteps=1000;
% pick some kind of time step?
dt=0.01;
% time array
t=zeros(1,nsteps);

x=zeros(nparticles,nsteps);
y=zeros(size(x));
z=zeros(size(x));
vx=zeros(size(x));
vy=zeros(size(x));
vz=zeros(size(x));
% the gravitational force has to be an array with nparticles elements?
gx=zeros(nparticles,1);
gy=zeros(nparticles,1);
gz=zeros(nparticles,1);
g=zeros(1,nsteps);
rmin=zeros(1,nsteps);

% spin array? keeps track of the inherent angular momentum of each particle
% For a start, I will initialize all particles with spin of zero, and then 
% only update inherent spin if two particles merge together, adding the angular
% momentum of each particle in the two-particle system.
s=zeros(nparticles,1);
% this array keeps track of the radius of each particle.
% In this simulation set up, the initial radius is a function of the initial particle density
% alpha_d is some dimensionless factor, less than one, that is the ratio of the
% particle radius to the initial inter-particle spacing 
% try 1/100 as a start
alpha_d=1/100;
% the characteristic length scale of the simulation volume, whether slab or 
% spherical.
lsim=10;
% the initial radius for each particle
rpart=alpha_d*lsim/nparticles.^(1/3);
% intial particle radius for all particles starts off the same
darr=ones(nparticles,1)*rpart; 
% mass array? Keeps track of different particles, maybe this could be unnecessary?
m=ones(nparticles,1).*(4/3*pi*rho*darr.^3);

[x(:,1),y(:,1),z(:,1),vx(:,1),vy(:,1),vz(:,1)]= ...
  particle_init(x(:,1),y(:,1),z(:,1),vx(:,1),vz(:,1),vy(:,1),lsim);
%rmin(1)=

% simple field calculator
[gx,gy,gz]=simple_g_field(x(:,1),y(:,1),z(:,1),m);
g(1)=max(sqrt(gx.^2+gy.^2+gz.^2));


% some tracer particles:
indices=randperm(nparticles);
%indices=indices(1:50);

%x(:,1)=[-1; 1];
%y(:,1)=[0; 0];
%z(:,1)=[0; 0];
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Go backward one half timestep, so that we have velocities at the t_(-1/2)
% time step.
% The particle velocities will be replaced with their t_(-1/2) values.
% I am pretty sure that xtemp, etc. should be equal to x(:,1), and so on
[xtemp,ytemp,ztemp,vx(:,1),vy(:,1),vz(:,1)]= ... 
      pusher(-dt/2,x(:,1),y(:,1),z(:,1),vx(:,1),vy(:,1),vz(:,1),gx,gy,gz,m,darr);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for k=2:nsteps
	
   % % display progress.
    if ~mod(k-2,nsteps/20)
		disp([num2str((k-2)/nsteps*100,'%.2f') '%'])
    end
    
    % where do I check to see if two particles have gotten close enough to merge?
    % probably in the pusher, i.e., at the half step, not the full step
    
    [x(:,k),y(:,k),z(:,k),vx(:,k),vy(:,k),vz(:,k)]= ... 
      pusher(dt,x(:,k-1),y(:,k-1),z(:,k-1),vx(:,k-1),vy(:,k-1),vz(:,k-1),gx,gy,gz,m,darr);
      
    %x(:,k)=xtemp;
    %y(:,k)=ytemp;
    %z(:,k)=ztemp;
    %vx(:,k)=vxtemp;
    %vy(:,k)=vytemp;
    %vz(:,k)=vztemp;
    
    t(k)=t(k-1)+dt;
    % compute new fields at the new particle positions
    
    % simple field calculator
    [gx,gy,gz]=simple_g_field(x(:,k),y(:,k),z(:,k),m);
    %g(k)=max(sqrt(gx.^2+gy.^2+gz.^2));
    %gtemp(k)=gx;
    % keep track of minimum distance between two particles?
    %rmin(k)=min(sqrt(x(:,k).^2
    
    % 
    %disp(num2str(k))
    %figure(1);clf;
    %drawnow;
    %scatter3(x(indices,k),y(indices,k),z(indices,k));
end

% make an animation?
bound=10;
for index=1:5:nsteps
  drawnow;
  scatter3(x(:,index),y(:,index),z(:,index),2)
  xlim([-bound bound]);
  ylim([-bound bound]);
  zlim([-bound bound]);
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~