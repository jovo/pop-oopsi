function Z=spPY(P,Y,X,t)
%HMM observation function P(Y|X,t) 
%Z=PY(p,[],X,t) computes a sample Y from P(Y|X,t). If X is a NxM matrix of 
%M dimensional states, Z is a NxK vector of the samples.
%Z=PY(p,Y,X,t) computes P(Y|X,t). If X is a NxM matrix of N M-dimensional 
%states, Z is a Nx1 vector of the probabilities (Y should be a scalar or 
%nan to indicate missing observations. "p" is a structure of parameters. 
%Yuriy Mishchenko 2009 Columbia Un
%
%REQUIRED PARAMETERS/DEFAULTS
% P.alpha  = 1;     %fluorescence scale
% P.beta   = 0;     %fluorescence offset
% P.gamma  = 1e-4;  %fluorescence conversion factor
% P.zeta   = 4e-4;  %fluorescence background noise
% P.n      = 1;     %fluorescence saturation exponent
% P.k_d    = 100;   %fluorescence saturation offset


% **************  SAMPLE  *****************  %
if(isempty(Y))                    %DRAW FROM P(Y|X,t) 
  C=X(:,2);                       %[Ca]
  S=Hill_v1(P,C);
  s2=P.gamma*S+P.zeta;            %variance
  S=P.alpha*S+P.beta;             %mean  
  Z=S+sqrt(s2).*randn(size(S));   %observation
else                              %CALCULATE P(Y|X,t)  
  C=X(:,2);                       %[Ca]  
  
  S=Hill_v1(P,C);
  s2=P.gamma*S+P.zeta;            %variance, expected
  S=P.alpha*S+P.beta;             %mean, expected    
  dY=(Y-S).^2./(2*s2);            %difference  
  adj=min(dY);                    %adjust by this constant
  Z=exp(adj-dY)./sqrt(s2);
  
  if(size(X,2)>=3) Z=Z./X(:,3); end%rho probabilities
end


function F = Hill_v1(P,C)
% generalized Hill model
C(C<0) = 0; F = C.^P.n./(C.^P.n+P.k_d);