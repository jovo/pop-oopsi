function [Z p]=PG_GOOPSI(P,X,t)
%Grid function rho(X,t) [for GOOPSI:Vogelstein-Paninski setup]
%Z=PG(P,[],t) returns a sample of X from rho(X,t).
%Z=PG(p,X,t) computes rho(X,t). If X is a NxM matrix of N M-dimensional
%states, Z is Nx1 vector of such probabilities P_i=rho(X_i,t).
%[Z,P]=PG(p,[],t) returns a sample of X from rho(X,t) and P=rho(X,t).
%"p" is a structure of parameters. Yuriy Mishchenko 2009 Columbia Un
%
%REQUIRED PARAMETERS
% P.nf              %1xT cell array of KxM arrays of M-dimensional particles.
% P.wf              %KxT array of particle weights.
% P.vf              %2xT array of kernel variances for [n C]'xT.
%OR
% P.alpha  = 1;     %fluorescence scale
% P.beta   = 0;     %fluorescence offset
% P.gamma  = 1e-4;  %fluorescence conversion factor
% P.zeta   = 4e-4;  %fluorescence background noise
% P.n      = 1;     %fluorescence saturation exponent
% P.k_d    = 100;   %fluorescence saturation offset
% P.F               %1xT array of observed fluorescence
% P.minVar = 1;     %lower limit on variance of C

if(~isfield(P,'nf')) %if swarm is not defined, use PRIOR SAMPLER
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
else %if swarm is defined, use KERNEL DENSITY ESTIMATION
  if(isempty(X))                  %DRAW SAMPLE FROM rho(X,t)
    %THIS IS THE ONLY PART THAT SHOULD RUN WITH MCMC SAMPLER
    ph=P.wf(:,t);                 %weights of different gaussians
    phc=cumsum(ph)/sum(ph);       %normalization
    r=rand; i=find(r<phc,1);      %draw the kernel from the mixture

    n=P.nf{t}(i,1); C=P.nf{t}(i,2);%unwrap the kernel from the mixture

    n=(rand>P.vf(1,t))==(n>0);    %draw rndized spike-state from kernel, TRICK HERE
    C=C+P.vf(2,t)*randn;          %draw rndized Ca-state from kernel
    Z=[n C];

    %EVALUATE rho(X,t)=P(n)*P(C)
    if(nargout==1) p=NaN; return; end    
    p1=sum(ph(P.nf{t}(:,1)>0));   %prob to pull n==1 from diff kernels
    %P(n):
    %P(1)=P(base=1)*P{stay}+P(base=0)*P{switch}
    %P(0)=P(base=1)*P{switch}+P(base=0)*P{stay}
    if(n>0) p=(1-P.vf(1,t))*p1+P.vf(1,t)*(1-p1); end
    if(n==0) p=P.vf(1,t)*p1+(1-P.vf(1,t))*(1-p1); end
    %P(C):
    %prob to pull C from diff Gauss kernels
    dC=(C-P.nf{t}(:,2)).^2;       
    s2=P.vf(2,t)^2;
    p1=exp(-dC/(2*s2))./sqrt(2*pi*s2);
    p=p*sum(phc.*p1);
  else                            %CALCULATE rho(X,t)
    Z=zeros(size(X,1),1);         %rho(X,t)=P(n)*P(C)
    ph=P.wf(:,t);                 %weights of different gaussians
    phc=cumsum(ph)/sum(ph);       %normalization
    for k=1:length(Z)
      n=X(k,1); C=X(k,2);         %unwrap state
      p1=sum(ph(P.nf{t}(:,1)>0)); %prob to pull n=1 from diff kernels
      %P(n):
      %P(1)=P(base=1)*P{stay}+P(base=0)*P{switch}
      %P(0)=P(base=1)*P{switch}+P(base=0)*P{stay}
      if(n>0) p=(1-P.vf(1,t))*p1+P.vf(1,t)*(1-p1); end
      if(n==0) p=P.vf(1,t)*p1+(1-P.vf(1,t))*(1-p1); end
      %P(C):
      %prob to pull C from diff Gauss kernels
      dC=(C-P.nf{t}(:,2)).^2;     
      s2=P.v(2,t)^2;
      p1=exp(-dC/(2*s2))./sqrt(2*pi*s2);
      Z(k)=p*sum(phc.*p1);
    end
    p=NaN;
  end
end


function C = AHill_v1(P,F)
% generalized hill model, inverse
C = (F./(1-F)).^(1/P.n)*P.k_d; C(C<0)  = 0;
