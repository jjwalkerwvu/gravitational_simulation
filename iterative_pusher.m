% % iterative_pusher.m
% %
% % iterative_pusher.m was written to treat drag forces on a dust grain. 
% % The Boris algorithm is incapable of treating drag terms or forces that 
% % are velocity^2 dependent. 
% % 7/18/2013
% %
% % currently have problems if the inputs vx, vy are both zero; need to 
% % fix this issue. Sept 2013.

% below is the version of inputs that I intend to use. Also, got rid of 
% global variables in here. (sept 2013)
function [x,y,vx,vy,w]=iterative_pusher(dtNwt,a,rho_d,q,x0,y0,vx0,vy0,...
    species,Ex,Ey,B,gx,gy,ne,ni,n_neut,vnx,vny,vex,vey,vix,viy,Te,Ti,...
    lambda_D,ch_model)
% Use ch_model flag to determine what description of drag force you want to
% use; this still needs to be built in to the code.

% % I've decided to get rid of global vars; they are commented if you feel 
% % like using them again.
%global qe;
%global me;
%global mp;
%global eps0;
qe=1.6e-19;
me=9.1e-31;
mp=1.67e-27;
eps0=8.854e-12;
mi=species*mp;

eta=ne/ni;
Tau=Te/Ti;
mr=me/mi;
vthe=sqrt(2*qe*Te/me);	% local electron thermal speed, m/s
vthi=sqrt(mr/Tau)*vthe;	% local ion (proton) thermal speed, m/s


m_neut=mi;
Tn=Ti;
C=4*pi*eps0*a;
%w=sqrt(vx0.^2+vy0.^2);
% % compute grain speed relative to an electron flow:
we=sqrt((vx0-vex).^2+(vy0-vey).^2);
% % compute grain speed relative to an ion flow:
% % velocity of the grain is relative to the ions:
% % e.g., if vx = 0, but the ions are streaming towards the grain in the 
% % -x direction with velocity vix, then it is equivalent to the grain 
% % moving at a velocity +vix in the +x-direction.
wi=sqrt((vx0-vix).^2+(vy0-viy).^2);
% % make a w-vector; the first element is the grain speed relative to 
% % electron flow, the second is the grain speed relative to ion flow.
w=[we wi];
% % compute some relevant parameters:
cn=sqrt(8*qe*Tn/pi/m_neut);
ci=sqrt(8*qe*Ti/pi/mi);
B0=abs(B);
md=4/3*pi*rho_d*a.^3;


% % Figure out whether to use Epstein drag force, or if it should go as 
% % velocity^2.
if wi<cn
	delta=1.26;
	% % nu_dn is the dust-neutral collision frequency; a result from the 
    % % Epstein drag
	nu_dn=delta*n_neut*cn*m_neut/a/rho_d;
	%dtNwt
	%pause
	% % line below is an older attempt at dust-neutral collision frequency; 
	% % vestigal and does not work presently.
	%nu_dn=4/3*pi*delta*n_neut*cn*m_neut*a.^2;
	beta_n=0;
	%disp('Epstein Drag')
	
	
	
else	
	% % nu_dn is the dust-neutral collision frequency; a result from the 
    % % Epstein drag
	nu_dn=0;
	% % think about if this can be rewritten (august 2013)
	beta_n=pi*a.^2*m_neut*n_neut*cn/md;
end
%nu_dn=0;
% % Set the ux=vx(t-1/2*dt) and uy=vy(t-1/2*dt), which are inputs
ux=vx0;
uy=vy0;	
vx=ux;
vy=uy;
% % Define some useful constants, which are coefficients for drag terms
% % (N.B. - They are not all in the same units!!!)
beta_ic=pi*ni*mi/md;
beta_io=4*pi*ni*mi/md;
% for testing
%beta_ic=0;beta_io=0;

% % Everything below here is needed in the main iterative loop.
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % max number of iterations, to prevent infinite loop
Nmax=1e3;
% % initialize the iteration counter to zero
n_iter=0;
% % Not sure how to initialize the error
err=1;
% % error tolerance; this should be scaled by some characteristic velocity.
tol=1e-10;

while err>tol && n_iter<=Nmax
	% % some useful definitions:
	vxdrifti=(vx+ux)/2-vix;
	vydrifti=(vy+uy)/2-viy;
	vxdriftn=(vx+ux)/2-vnx;
	vydriftn=(vy+uy)/2-vny;
	vs=sqrt(8*qe*Ti/pi/mi+vxdrifti^2+vydrifti^2);
	vn=sqrt(vxdriftn.^2+vydriftn.^2);
	bc=sqrt(1-2*qe*q/C/mi/(vs.^2))*a;
    %b_90=qe*abs(q)/(4*pi*eps0*mi*vs.^2);
	b_90=qe*q/(4*pi*eps0*mi*vs.^2);     % collision paramter for 90 degree 
                                        % collisions
	gamma = 1/2*log((lambda_D^2+b_90^2)/(bc.^2+b_90^2));	
    % see Grabbe-Merlino or other Dusty plasma texts for more information.


	% % Elements of the Jacobian matrix, which are derivatives of f1 and f2 
    % % with respect to vx, vy in matrix form:
	% %	| (df1/dvx) 	(df1/dvy) 	|(vx^[k+1] - vx^[k])	=-(f1)
	% %	| (df2/dvx) 	(df2/dvy) 	|(vy^[k+1] - vy^[k])	=-(f2)
	b=zeros(2,1);
	ut=zeros(2,1);
	lt=zeros(2,1);
	f=zeros(2,1);
	
	% % make sure to check if vn=0!
	if vn==0
		b(1)=-1+.5*dtNwt*(-nu_dn-beta_n*(vn+(vxdriftn)^2/vn)-...
            beta_ic*vs*bc.^2-beta_io*vs*b_90.^2*gamma...
		-beta_ic*(vxdrifti).^2*(bc.^2/vs+2*vs*a.^2*(2*qe*q/C/mi)/vs.^4)...
		-beta_io*(vxdrifti).^2*(b_90.^2*gamma/vs-4*b_90.^2*gamma/vs...
		+2*b_90.^2*exp(-gamma)*(b_90.^2*(a.^2*2*qe*q/C/mi/vs.^2-bc.^2)...
		+lambda_D.^2*(b_90.^2-2*a.^2*qe*q/C/mi/vs.^2))/vs/(bc.^2+b_90.^2).^2));
		
		b(2)=0;
		ut(1)=0;
		lt(2)=0;
		disp('vn = 0')
	% % Proceed as normal if vn ~= 0.	
	else
		% % df1/dvx:
		b(1)=-1+.5*dtNwt*(-nu_dn-beta_n*(vn+(vxdriftn)^2/vn)-...
            beta_ic*vs*bc.^2-beta_io*vs*b_90.^2*gamma...
		-beta_ic*(vxdrifti).^2*(bc.^2/vs+2*vs*a.^2*(2*qe*q/C/mi)/vs.^4)...
		-beta_io*(vxdrifti).^2*(b_90.^2*gamma/vs-4*b_90.^2*gamma/vs...
		+2*b_90.^2*exp(-gamma)*(b_90.^2*(a.^2*2*qe*q/C/mi/vs.^2-bc.^2)...
		+lambda_D.^2*(b_90.^2-2*a.^2*qe*q/C/mi/vs.^2))/vs/(bc.^2+b_90.^2).^2));
		
		% % df2/dvy:
		b(2)=-1+.5*dtNwt*(-nu_dn-beta_n*(vn+(vydriftn)^2/vn)-...
            beta_ic*vs*bc.^2-beta_io*vs*b_90.^2*gamma...
		-beta_ic*(vydrifti).^2*(bc.^2/vs+2*vs*a.^2*(2*qe*q/C/mi)/vs.^4)...
		-beta_io*(vydrifti).^2*(b_90.^2*gamma/vs-4*b_90.^2*gamma/vs...
		+2*b_90.^2*exp(-gamma)*(b_90.^2*(a.^2*2*qe*q/C/mi/vs.^2-bc.^2)...
		+lambda_D.^2*(b_90.^2-2*a.^2*qe*q/C/mi/vs.^2))/vs/(bc.^2+b_90.^2).^2));
		
		% % "upper" triangular portion of the Jacobian matrix, or df1/dvy:
		% % (In tridiagonal terms, this is the "C" vector term)
		%df1_dvy=ut(1)
		ut(1)=.5*dtNwt*(q*B0/md-nu_dn-beta_n*(vn+(vydriftn)^2/vn)-...
            beta_ic*vs*bc.^2-beta_io*vs*b_90.^2*gamma...
		-beta_ic*(vydrifti).^2*(bc.^2/vs+2*vs*a.^2*(2*qe*q/C/mi)/vs.^4)...
		-beta_io*(vydrifti).^2*(b_90.^2*gamma/vs-4*b_90.^2*gamma/vs...
		+2*b_90.^2*exp(-gamma)*(b_90.^2*(a.^2*2*qe*q/C/mi/vs.^2-bc.^2)...
		+lambda_D.^2*(b_90.^2-2*a.^2*qe*q/C/mi/vs.^2))/vs/(bc.^2+b_90.^2).^2));

			
	
		% % "lower" triangular portion of the Jacobian matrix or df2/dvx:
		% % (In tridiagonal terms, this is the "A" vector term)
		%df2_dvx=lt(2)
		lt(2)=.5*dtNwt*(-q*B0/md-nu_dn-beta_n*(vn+(vxdriftn)^2/vn)-...
            beta_ic*vs*bc.^2-beta_io*vs*b_90.^2*gamma...
		-beta_ic*(vxdrifti).^2*(bc.^2/vs+2*vs*a.^2*(2*qe*q/C/mi)/vs.^4)...
		-beta_io*(vxdrifti).^2*(b_90.^2*gamma/vs-4*b_90.^2*gamma/vs...
		+2*b_90.^2*exp(-gamma)*(b_90.^2*(a.^2*2*qe*q/C/mi/vs.^2-bc.^2)...
		+lambda_D.^2*(b_90.^2-2*a.^2*qe*q/C/mi/vs.^2))/vs/(bc.^2+b_90.^2).^2));
	
		
	end
	
		
	% % The functions of vx, vy:
	f(1)=-vx+ux+dtNwt*(gx+q*Ex/md+q*B0*(vy+uy)/2/md-nu_dn*...
        (vxdriftn)-beta_n*(vxdriftn)*vn...
		-beta_ic*vs*(vxdrifti)*bc.^2-beta_io*vs*(vxdrifti)*b_90.^2*gamma);
	f(2)=-vy+uy+dtNwt*(gy+q*Ey/md-q*B0*(vx+ux)/2/md-nu_dn*(vydriftn)-...
        beta_n*(vydriftn)*vn...
		-beta_ic*vs*(vydrifti)*bc.^2-beta_io*vs*(vydrifti)*b_90.^2*gamma);
	
	% % Need the negative of f:
	f=-f;
	
	% % the "error" vector can be computed using a tridiagonal solver.
	error=trisolver(lt,b,ut,f);
	
	% % figure out the maximum error in the "error vector", if this is less 
	% % than the tolerance then break out of the loop. Use the "infinity" 
	% % norm, or maximum value or the error vector.
	err=max(abs(error));
	vx=vx+error(1);
	vy=vy+error(2);
	% % put this line below to prevent an infinite loop
	n_iter=n_iter+1;
	
	%q/1.6e-19
	%pause
end
% % % velocity of the grain relative to the ions:
% % % e.g., if vx = 0, but the ions are streaming towards the grain in the 
% % % -x direction with velocity vix, then it is equivalent to the grain 
% % % moving at a velocity +vix in the +x-direction.
% %w=sqrt((vx-vix).^2+(vy-viy).^2);
% % % compute grain speed relative to an electron flow:
% we=sqrt((vx-vex).^2+(vy-vey).^2);
% % % compute grain speed relative to an ion flow:
% % % velocity of the grain relative to the ions:
% % % e.g., if vx = 0, but the ions are streaming towards the grain in the 
% % % -x direction with velocity vix, then it is equivalent to the grain 
% % % moving at a velocity +vix in the +x-direction.
% wi=sqrt((vx-vix).^2+(vy-viy).^2);
% % % make a w-vector; the first element is the grain speed relative to 
% % % electron flow, the second is the grain speed relative to ion flow.
% w=[we wi];

% In a sheath, we assume that the electrons are not flowing. We will allow
% for electron drifts in the phi direction, which requires that vey is set
% up in profiles properly so that vey = ve_phi.
we = abs(vey);
% % compute grain speed relative to an ion flow:
% % velocity of the grain relative to the ions:
% % e.g., if vx = 0, but the ions are streaming towards the grain in the 
% % -x direction with velocity vix, then it is equivalent to the grain 
% % moving at a velocity +vix in the +x-direction.
% WE ADD THE TERM (vex).^2 IN ORDER TO GET THE RELATIVE DRIFT ALONG THE 
% Z-DIRECTION. THE DUST GRAIN, BY VIRTUE OF BEING LEVITATED IN THE PLANAR
% SHEATH, HAS NO VELOCITY ALONG THE Z DIRECTION.
wi=sqrt((vx-vix).^2+(vy-viy).^2 + (vex).^2);
% % make a w-vector; the first element is the grain speed relative to 
% % electron flow, the second is the grain speed relative to ion flow.
w=[we wi];


%%~~~~#2
% % second main step in this iterative method: calculate positions at the 
% % full timestep based on the velocities calculated at the half timestep.  
% % (The positions will be half a timestep ahead of the velocities). This 
% % is done after vx and vy have been found through an iterative process.


x=dtNwt*vx+x0;
y=dtNwt*vy+y0;
	
end