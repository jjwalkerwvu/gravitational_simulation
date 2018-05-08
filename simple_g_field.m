% simple_g_field.m

% My lowest order attempt at making a gravitational field calculator; just sum 
% up the forces on each particle from every other particle

function [gx,gy,gz]=simple_g_field(x,y,z,marr);

%G=6.7e-11;
G=1;
m=1;


for index=1:length(x);
 
  % the (x~=x(index)) indexing is to prevent us from finding the distance 
  % between the "indexed" particle and itself, which we already know is zero
  mvect=marr;
  % discard the value in the mass matrix corresponding to the current particle,
  % upon which we are calculating all the forces from other particles
  mvect(index)=[];
  vectx=(x(x~=x(index))-x(index));
  vecty=(y(y~=y(index))-y(index));
  vectz=(z(z~=z(index))-z(index));
  dist=sqrt((x(x~=x(index))-x(index)).^2+(y(y~=y(index))-y(index)).^2+(z(z~=z(index))-z(index)).^2);
  % find particles that are "touching?"
  % remove possible divisions by zero
  
  % Do I still need these in case two particles are somehow touching? 
  % I think not, because whenever they get close enough, they should be merged
  %remove_array=find(dist==0);
  %vectx(remove_array)=[];
  %vecty(remove_array)=[];
  %vectz(remove_array)=[];
  %dist(remove_array)=[];
  % need to somehow make particles capture each other when they get too close
  %if min(dist)<5e-2;
  %  disp("too close");
  %end
  %if length(nnz(dist))<length(dist)
  %  disp(num2str(gx(index))
  %  5+5;
  %end
  % need negative signs in front of all these vector quantities, because gravity
  % is an attractive force
  gx(index)=G*marr(index)*sum(mvect.*vectx./dist.^3);
  gy(index)=G*marr(index)*sum(mvect.*vecty./dist.^3);
  gz(index)=G*marr(index)*sum(mvect.*vectz./dist.^3);
end

% I think we need to transpose gx, gy, gz arrays so that they are in the correct
% format.

gx=gx';
gy=gy';
gz=gz';