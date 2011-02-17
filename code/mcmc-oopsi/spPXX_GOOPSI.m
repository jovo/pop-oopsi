function Z=spPXX_GOOPSI(P,X1,X0,t)
%HMM transition function P(X1|X0,t) [for GOOPSI:Vogelstein-Paninski setup] 
%Z=PXX(p,[]) computes a 1xM-dimensional sample from P(X,t==0).
%Z=PXX(p,X1) computes the probability P(X1,t==0). If X1 is NxM matrix of N 
%M-dimensional states, Z is a Nx1 vector of such probabilities.
%Z=PXX(p,[],X0,t) computes a 1xM-dimensional sample from P(X1|X0,t). 
%Z=PXX(p,X1,X0,t) computes P(X1|X0,t). If X1 is NxM matrix and X0 is KxM 
%matrix of N and K M-dimensional states, respecitvely, Z is a NxK matrix
%of transition probabilities between these states P(i,j)=P(X1(i)|X0(j),t). 
%"p" is a structure of parameters. Yuriy Mishchenko 2009 Columbia Un
%
%REQUIRED PARAMETERS/DEFAULTS
% P.dt      = 0.02;    %time bin
% P.k       = 1;       %spike generation probabilities, may be 1xT array
% P.A       = 50;      %Ca jump 
% P.sigma_c = 15;      %Ca noise
% P.tau_c   = 0.2;     %Ca reporter time constant
% P.C_0     = 50;      %Ca background
%[optional]
% P.sig2    = P.sigma_c^2*P.dt;%variance adjusted for dt
% P.rate    = exp(-P.k*P.dt);%spiking rate, may be 1xT array
if(nargin<4) t=1; end

if(isfield(P,'sig2')) s2=P.sig2; else s2=P.sigma_c^2*P.dt; end%precompute
if(isfield(P,'rate'))                                         %precompute
  if(length(P.rate)>=t) p=P.rate(t); else p=P.rate(end); end
else
  if(length(P.k)>=t) p=exp(-P.k(t)*P.dt); else p=exp(-P.k(end)*P.dt); end
end

if(nargin<3 && isempty(X1))             %DRAW FROM P(X,t==0)
  n=rand>p;                             %spike
  C=P.C_0+sqrt(s2)*randn;               %base Ca
  Z=[n C];
elseif(nargin<3)                        %CALCULATE P(X,t==0)
  %THIS IS ONE PART THAT SHOULD RUN WITH MCMC SAMPLER    
  n=X1(:,1); C=X1(:,2);                 %unwrap state
  dC=(C-P.C_0).^2/(2*s2);               %expected Ca-level
  adj=min(dC(:));
  Z=exp(adj-dC);                        %base Ca prob
  Z(n>0) =Z(n>0)*(1-p);                 %spike prob
  Z(n==0)=Z(n==0)*p;
elseif(isempty(X1))                     %DRAW FROM P(X1|X0,t)
  n=X0(1); C=X0(2);                     %unwrap previous point
  n=rand>p;                             %new spike    
  C=C-P.dt/P.tau_c*(C-P.C_0)+P.A*n+sqrt(s2)*randn;%new Ca
  Z=[n C];
else                                    %CALCULATE P(X1|X0,t)
  %THIS IS THE PART THAT SHOULD RUN WITH MCMC SAMPLER  
  n0=X0(:,1); C0=X0(:,2); M=length(C0);
  n1=X1(:,1); C1=X1(:,2); N=length(C1);

  warning off all
  C0=C0-P.dt/P.tau_c*(C0-P.C_0);        %expected [Ca]-left
  C1=C1-P.A*n1;                         %expected [Ca]-right
  A=repmat(C1(:),[1 M]);
  B=repmat(C0(:),[1 N]);
  A=(A-B').^2/(2*s2);                   %deviations
  %PRECONDITIONER
  if(isfield(P,'prc')&&~isempty(P.prc)) PRC=log(P.prc); else PRC=0; end
  A=A-PRC;   
  adj=min(A(:));                        %adjust Pr by this constant
  Z=exp(adj-A);                         %[Ca]-[Ca] transition prob  
  Z(isnan(Z))=0;
  warning on all

  Z(n1>0,:) =Z(n1>0,:)*(1-p);           %n1 generation prob
  Z(n1==0,:)=Z(n1==0,:)*p;
end