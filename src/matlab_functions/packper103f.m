function [x y D Gn nu Ek Ep Ct Dt nt tc]=packper103f(Lx,Ly,Dn,Gam,seed,ispack,plotit,tol,ic,fx,fy)
rng(seed)
% Example N = 10; Dn_list = 1+0.4*(rand([1,N])>0.5); packper103f(1,1,Dn_list,0,1)
% Molecular Dynamics Simulator
% <md.m> Mark D. Shattuck 7/22/2010

% revision history:
% 07/22/2010 Mark D. Shattuck <mds> md.m
%            MD Demo for HandsOn 2010
%            000 mds Initial conditions and visualization
% 07/24/2010 001 mds Add Euler
% 07/24/2010 002 mds Add Interaction detection and Force Law
% 07/24/2010 003 mds Add Velocity Verlet
% 07/24/2010 004 mds Add Nplotskip
% 08/22/2010 006 mds Add saving
% 10/11/2010 007 mds Add add fast force
% 11/18/2010 000 mds pack.m Packing finder redo of fastpack
% 04/18/2012 100 mds packper.m periodic packing branch 1
% 05/04/2012 100 mds function version
% 02/25/2013 101 mds new hunt and stop from polypack
% 11/14/2014 103 mds remove floater calc for speed

%% Experimental Parmeters


N=length(Dn);    % number particles
Ds=min(Dn);      % smallest particle
%Dl=max(Dn);      % largest particle
Gn=Dn/Ds;        % particle ratios
K=100;   % spring constant for harmonic force law
B=1;     % damping
P = 0;
% File naming - just numbers to fit into the simulation_2D.m function for testing.
P_target = 0; % big enough to know that this is NOT a regular packing
W_factor = 0;
filename = ['in/2D_Repeating/2D_N' num2str(N) '_P' num2str(P_target) '_Width' num2str(W_factor) '_Seed' num2str(seed) '.mat'];
% %%%%%

if(~exist('tol','var')|| isempty(tol))
  tol=1e-6;
end

Ptol=tol;    % potential energy per particle<Ptol in final state
Pthresh=1e-5; % potential per particle<Pthresh before growth
Pstep=1e-4;   % potential per particle for add in growth

if(~exist('ispack','var')|| isempty(ispack))
  ispack=true;
end

if(~exist('ic','var') || isempty(ic))
  ic='Rand';
end

%% Display Parameters
if(~exist('plotit','var') || isempty(plotit))
  plotit=true;
end
Nplotskip=200;  % number of timesteps to skip before plotting

%% Simulation Parmeters
dt=.5e-1;        % time step
Nt=200000;       % max time per pack

%% Calculated Parameters
Gnm=mdiff(Gn,-Gn)/2;  % min separation matrix ratio Dnm=Ds*Gnm

%% Initial Conditions
switch ic
  case 'Rand'
    x=Lx*rand(1,N);
    y=Ly*rand(1,N);
  case 'TimeTest'
    load('dilute100','x','y','D','N','Lx','Ly');
    x=x+randn(RandStream('mt19937ar','seed',100),1,N)./Ds/10; %#ok<NODEF>
    y=y+randn(RandStream('mt19937ar','seed',200),1,N)./Ds/10; %#ok<NODEF>
  case 'Fix'
    x=fx;
    y=fy;
end

x0=x; %#ok<NASGU>
y0=y; %#ok<NASGU>

%% find initial D
dx=repmat(x,N,1);
dx=dx-dx';
dx=dx-round(dx/Lx)*Lx;  % Periodic x
dy=repmat(y,N,1);
dy=dy-dy';
dy=dy-round(dy/Ly)*Ly;  % Periodic x
dnm=sqrt(dx.^2+dy.^2); % Distance between all particles
dnm(1:N+1:end)=inf;    % set diagonal to D
D=min(dnm(:)./Gnm(:));

Dn=D*Gn;        % Diameter list

vx=zeros(1,N);
vy=zeros(1,N);

ax_old=0*x;
ay_old=0*y;


%% Save Variables
Ek=zeros(Nt,1);    % Kinetic Energy
Ep=zeros(Nt,1);    % particle-particle potential
Ct=zeros(Nt,1);    % particle-particle potential
%Cft=zeros(Nt,1);    % particle-particle potential
Dt=zeros(Nt,1);    % particle-particle potential
%Nft=zeros(Nt,1);    % particle-particle potential

%% Setup Plotting
if(plotit)
  clf;
  h=zeros(1,4*N);
  for np=1:N
    h(np)=rectangle('Position',[x(np)-.5*Dn(np) y(np)-.5*Dn(np) Dn(np) Dn(np)],'Curvature',[1 1],'edgecolor','b');
    h(np+N)=rectangle('Position',[2 2 Dn(np) Dn(np)],'Curvature',[1 1],'edgecolor','b');
    h(np+2*N)=rectangle('Position',[2 2 Dn(np) Dn(np)],'Curvature',[1 1],'edgecolor','b');
    h(np+3*N)=rectangle('Position',[2 2 Dn(np) Dn(np)],'Curvature',[1 1],'edgecolor','b');
  end
  Np=N;
  
  % Left wall
  ii=find(x<Dn/2);
  for nn=1:length(ii)
    np=ii(nn);
    set(h(nn+Np),'Position',[x(np)+Lx-.5*Dn(np) y(np)-.5*Dn(np) Dn(np) Dn(np)]);
  end
  Np=Np+length(ii);
  
  % Bottom wall
  ii=find(y<Dn/2);
  for nn=1:length(ii)
    np=ii(nn);
    set(h(nn+Np),'Position',[x(np)-.5*Dn(np) y(np)-.5*Dn(np)+Ly Dn(np) Dn(np)]);
  end
  Np=Np+length(ii);
  
  %Right wall
  ii=find(x>Lx-Dn/2);
  for nn=1:length(ii)
    np=ii(nn);
    set(h(nn+Np),'Position',[x(np)-.5*Dn(np)-Lx y(np)-.5*Dn(np) Dn(np) Dn(np)]);
  end
  Np=Np+length(ii);
  
  % Top wall
  ii=find(y>Ly-Dn/2);
  for nn=1:length(ii)
    np=ii(nn);
    set(h(nn+Np),'Position',[x(np)-.5*Dn(np) y(np)-.5*Dn(np)-Ly Dn(np) Dn(np)]);
  end
  Np=Np+length(ii);
  
  % Corner
  dx=x-round(x/Lx)*Lx;
  dy=y-round(y/Ly)*Ly;
  ii=find(sqrt(dx.^2+dy.^2)<Dn/2);
  for nn=1:length(ii)
    np=ii(nn);
    px=-sign(x(np)-Lx/2)*Lx;
    py=-sign(y(np)-Ly/2)*Ly;
    set(h(nn+Np),'Position',[x(np)-.5*Dn(np)+px y(np)-.5*Dn(np)+py Dn(np) Dn(np)]);
  end
  Np=Np+length(ii);
  
  for nn=1:4*N-Np
    set(h(Np+nn),'Position',[2 2 Dn(np) Dn(np)],'Curvature',[1 1],'edgecolor','b');
  end
  
  axis('equal');
  axis([0 Lx 0 Ly]);
  %pause;
end

%% Main Loop
Dh=-1;
Dl=-1;
gam=Gam;
%Cnm=zeros(N,N);
Cn=zeros(1,N);
for nt=1:Nt
  % plot particles
  if(plotit && (rem(nt-1,Nplotskip)==0))
    for np=1:N
      set(h(np),'Position',[x(np)-.5*Dn(np) y(np)-.5*Dn(np) Dn(np) Dn(np)]);
    end
    Np=N;
    
    % Left wall
    ii=find(x<Dn/2);
    for nn=1:length(ii)
      np=ii(nn);
      set(h(nn+Np),'Position',[x(np)+Lx-.5*Dn(np) y(np)-.5*Dn(np) Dn(np) Dn(np)]);
    end
    Np=Np+length(ii);
    
    % Bottom wall
    ii=find(y<Dn/2);
    for nn=1:length(ii)
      np=ii(nn);
      set(h(nn+Np),'Position',[x(np)-.5*Dn(np) y(np)-.5*Dn(np)+Ly Dn(np) Dn(np)]);
    end
    Np=Np+length(ii);
    
    %Right wall
    ii=find(x>Lx-Dn/2);
    for nn=1:length(ii)
      np=ii(nn);
      set(h(nn+Np),'Position',[x(np)-.5*Dn(np)-Lx y(np)-.5*Dn(np) Dn(np) Dn(np)]);
    end
    Np=Np+length(ii);
    
    % Top wall
    ii=find(y>Ly-Dn/2);
    for nn=1:length(ii)
      np=ii(nn);
      set(h(nn+Np),'Position',[x(np)-.5*Dn(np) y(np)-.5*Dn(np)-Ly Dn(np) Dn(np)]);
    end
    Np=Np+length(ii);
    
    % Corner
    dx=x-round(x/Lx)*Lx;
    dy=y-round(y/Ly)*Ly;
    ii=find(sqrt(dx.^2+dy.^2)<Dn/2);
    for nn=1:length(ii)
      np=ii(nn);
      px=-sign(x(np)-Lx/2)*Lx;
      py=-sign(y(np)-Ly/2)*Ly;
      set(h(nn+Np),'Position',[x(np)-.5*Dn(np)+px y(np)-.5*Dn(np)+py Dn(np) Dn(np)]);
    end
    Np=Np+length(ii);
    
    for nn=1:4*N-Np
      set(h(nn+Np),'Position',[2 2 Dn(np) Dn(np)],'Curvature',[1 1],'edgecolor','b');
    end
    title(num2str(pi/(4*Lx*Ly)*sum(Dn.^2)))
    drawnow;
  end
  
  Ek(nt)=sum((vx.^2+vy.^2))/2;
  
  x=x+vx*dt+ax_old.*dt.^2/2;  % first step in Verlet integration
  y=y+vy*dt+ay_old.*dt.^2/2;
  
  % Interaction detector and Force Law
  Fx=zeros(1,N);
  Fy=zeros(1,N);
  
 % Cnm(:)=0;
  Cn(:)=0;
  C=0;   % number of contacts
  Gmin=9e9;
  for nn=1:N
    for mm=nn+1:N
      dy=y(mm)-y(nn);
      im=round(dy/Ly);
      dy=dy-im*Ly;  % Periodic x
      Dnm=(Dn(nn)+Dn(mm))/2;
      if(abs(dy)<Dnm)
        dx=x(mm)-x(nn);
        dx=dx-round(dx/Lx-im*gam)*Lx-im*gam*Lx;
        dnm=sqrt(dx.^2+dy.^2);
        Gmin=min(Gmin,dnm/Dnm);
        if(dnm<Dnm)
          C=C+1;
          Cn(nn)=Cn(nn)+1;
          Cn(mm)=Cn(mm)+1;
%          Cnm(nn,mm)=1;
%          Cnm(mm,nn)=1;
          F=-K*(Dnm/dnm-1);
          Ep(nt)=Ep(nt)+K/2*(Dnm-dnm).^2;
          Fx(nn)=Fx(nn)+F.*dx;  % particle-particle Force Law
          Fx(mm)=Fx(mm)-F.*dx;
          Fy(nn)=Fy(nn)+F.*dy;  % particle-particle Force Law
          Fy(mm)=Fy(mm)-F.*dy;
        end
      end
    end
  end
  
  Fx=Fx-B*vx;
  Fy=Fy-B*vy;
  
%   done=false;
%   Nflt=0;
%   while(~done)
%     Cn=sum(Cnm);
%     iflt=find(Cn<2);
%     Cnm(iflt,:)=0;
%     Cnm(:,iflt)=0;
%     if(numel(iflt)>Nflt)
%       Nflt=numel(iflt);
%     else
%       done=true;
%     end
%   end
%   
%   Nft(nt)=Nflt;
%   Cft(nt)=sum(Cnm(:))/2;

  Ct(nt)=C;
  Dt(nt)=D;
 if (Ep(nt)<1e-10); 
   dt=1;
 end
 if(Ep(nt)>1e-8)  
   dt=.5e-1;
 end
 if(Ep(nt)<Pthresh*N && Dh<0)
    if(~ispack)
      if(D==Ds)
        break;
      end
    end
    D=D+sqrt(Pstep/K);
    Dn=D*Gn;
    if(~ispack)
      D=min(D,Ds);
    end
    %Dnm=D*Gnm;
  elseif(Ep(nt)<Ptol*N/100 && Dh>0)
    Dl=D;
    D=(Dl+Dh)/2;
    Dn=D*Gn;
  elseif(Ek(nt)<1e-15*N && C>N && Ep(nt)>Ptol*N)
    if(Dh<0)
      Dh=D;
      tc=nt;
    end
    if(Dl>0)
      Dh=D;
      D=(Dl+Dh)/2;
    else
      D=D-sqrt(Pstep/K);
    end
    Dn=D*Gn;
  elseif(Ek(nt)<1e-25*N && C>N)
    break;
  end
   
  ax=Fx;
  ay=Fy;
  
  vx=vx+(ax_old+ax).*dt/2;  % second step in Verlet integration
  vy=vy+(ay_old+ay).*dt/2;

  vx(Cn==0)=0;  % zero floater velocities
  vy(Cn==0)=0;
  
  ax_old=ax;
  ay_old=ay;
end

nu=pi/(4*Lx*Ly)*sum(Dn.^2);

% save(sprintf('%s%02d',basename,n,'.mat'),'x','y','D','Lx','Ly','nu','B','N','K','G','Ns');
fprintf(1,' %f\n',nu);

% Ek0=Ek(nt);
% Ep0=Ep(nt);
% C0=Ct(nt);
% 

save(filename, 'x', 'y', 'Dn', 'Lx', 'Ly', 'K', 'P_target', 'P');