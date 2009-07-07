function [Z p]=PG(P,X,t)
%Grid function rho(X,t)
%Z=PG(P,[],t) returns a sample of X from rho(X,t).
%Z=PG(p,X,t) computes rho(X,t). If X is a NxM matrix of N M-dimensional
%states, Z is Nx1 vector of such probabilities P_i=rho(X_i,t).
%[Z,P]=PG(p,[],t) returns a sample of X from rho(X,t) and P=rho(X,t).
%"p" is a structure of parameters. Yuriy Mishchenko 2009 Columbia Un
%
%REQUIRED PARAMETERS
% P.alpha  = 1;     %fluorescence scale
% P.beta   = 0;     %fluorescence offset
% P.gamma  = 1e-4;  %fluorescence conversion factor
% P.zeta   = 4e-4;  %fluorescence background noise
% P.n      = 1;     %fluorescence saturation exponent
% P.k_d    = 100;   %fluorescence saturation offset
% P.F               %1xT array of observed fluorescence
% P.minVar = 1;     %lower limit on variance of C


% **************  SAMPLE  *****************  %
p=0.5;                        %spike draw is 1/2
if(isempty(X))                %DRAW SAMPLE FROM rho(X,t)
  n=rand>p;                   %random spike
  S=(P.F(t)-P.beta)/P.alpha;  %expected S
  C=AHill_v1(P,S);            %expected C
  s=sqrt(P.gamma*S+P.zeta);   %expected variance
  s=2*C/(P.n*S*(1-S))*s;      %variance/Jacobian,
  s=max(P.minVar,s);          %but not too small
  dC=s*randn;
  Z=[n C+dC];                 %new state

  if(nargout==1) p=NaN; return; end
  p1=exp(-dC^2/(2*s^2))/(sqrt(2*pi)*s);%rho(x)=rho(n)*rho(C)
  if(n>0) p=p1*(1-p); else p=p1*p; end
else                          %CALCULATE rho(X,t)
  S=(P.F(t)-P.beta)/P.alpha;  %expected S
  C=AHill_v1(P,S);            %expected Ca
  s=sqrt(P.gamma*S+P.zeta);   %expected variance
  s=2*C/(P.n*S*(1-S))*s;      %variance/Jacobian,
  s=max(P.minVar,s);          %but not too small
  C0=X(:,2);                  %proposed Ca
  n=X(:,1);                   %proposed spike
  Z=exp(-(C0-C).^2/(2*s.^2))./(sqrt(2*pi)*s);%rho(x)=rho(n)*RHO(C)
  z=repmat(p,size(n)); z(n>0)=1-p;%rho(x)=RHO(n)*rho(C)
  Z=Z.*z;
  p=NaN;
end

function C = AHill_v1(P,F)
% generalized hill model, inverse
C = (F./(1-F)).^(1/P.n)*P.k_d; C(C<0)  = 0;
