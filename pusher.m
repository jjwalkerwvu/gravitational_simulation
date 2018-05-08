% pusher.m

% this is a simple particle pusher for moving particles in a gravitational field

function [x,y,z,vx,vy,vz]=pusher(dt,x0,y0,z0,vx0,vy0,vz0,gx,gy,gz,marr,rpart)

% m is the mass matrix; keeps track of each particle's mass
% rpart is an array containing the radius of each particle

% dv=dt*gx/m+dt*gy/m+dt*gz/m
% v=v0+dv?

% replaced with a mass matrix
%m=1;

vx=dt*gx./marr+vx0;
vy=dt*gy./marr+vy0;
vz=dt*gz./marr+vz0;


x=dt*vx+x0;
y=dt*vy+y0;
z=dt*vz+z0;

% after we have the new positions, check to see if any particles are "touching"
% one another; the matrix rpart keeps track of particle diameters

% what if all of this goes into the force calculation subroutine instead?
% positions and momenta from previous time step would be known, and it makes 
% sense to calculate whether particles touch on the whole timestep, which is
% also when new positions are calculated.
for index=1:length(x);
  % only find distances to other particles, not the current particle indexed
  % in the current loop iteration.
  % so you will have nparticles-1 array elements
  dist=sqrt((x(x~=x(index))-x(index)).^2+(y(y~=y(index))-y(index)).^2+ ...
    (z(z~=z(index))-z(index)).^2);
  % overlap is initialized as a dummy array
  overlap=rpart;
  % remove the array element corresponding to the indexed particle; we know that
  % it would overlap with itself.
  overlap(index)=[];
  overlap=dist<(rpart(index)+overlap);
  % need to take the transpose, to get this dummy array into the same dimension
  % as overlap array
  dummy=(1:length(x))';
  % again, get rid of the indexed particle
  dummy(index)=[];
  % multiply the logical array by the array containing 1,2,...,npart array to
  % get all the particles that overlap the current indexed particle
  overlap=dummy.*overlap;
  overlap=overlap(find(overlap));
end



% possibly merge particles, angular momenta, translational momenta after the 
% position update?

% if the particles are touching, use the position and momentum from previous time
% step to obtain the inherent spin of the newly-merged particle, and use 
% conservation of momentum to find the momentum of the new particle?
% or do I need to recalculate forces?