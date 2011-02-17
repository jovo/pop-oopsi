function [Z p]=PG_GOOPSI(P,X,t,N)
%Grid function rho(X,t,N) [for GOOPSI:Vogelstein-Paninski setup]
%Z=PG(P,[],t,N) returns a sample of N X's from rho(X,t).
%Z=PG(p,X,t) computes rho(X,t). If X is a NxM matrix of N M-dimensional
%states, Z is Nx1 vector of such probabilities P_i=rho(X_i,t).
%[Z,P]=PG(p,[],t,N) returns a sample of N X's from rho(X,t) and P=rho(X,t).
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
if(nargin<4 || isempty(N)) N=1; end

if(~isfield(P,'nf')) %if swarm is not defined, use PRIOR SAMPLER
  p=0.5;                        %prob no-spike 
  if(isempty(X))                %DRAW SAMPLE FROM rho(X,t)
    n=rand(N,1)>p;              %random spike
    S=(P.F(t)-P.beta)/P.alpha;  %expected S
    C=AHill_v1(P,S);            %expected C
    s=sqrt(P.gamma*S+P.zeta);   %expected variance
    s=2*C/(P.n*S*(1-S))*s;      %variance/Jacobian,
    s=max(P.minVar,s);          %but not too small
    dC=s*randn(N,1);
    Z=[n C+dC];                 %new state
    
    if(nargout==1) return; end
    p1=exp(-dC.^2/(2*s^2));     %rho(x)=rho(n)*rho(C)
    p1(n>0) =p1(n>0)*(1-p); 
    p1(n==0)=p1(n==0)*p;
    p=p1;
  else                          %CALCULATE rho(X,t)
    S=(P.F(t)-P.beta)/P.alpha;  %expected S
    C=AHill_v1(P,S);            %expected Ca
    s=sqrt(P.gamma*S+P.zeta);   %expected variance
    s=2*C/(P.n*S*(1-S))*s;      %variance/Jacobian,
    s=max(P.minVar,s);          %but not too small
    dC=X(:,2)-C;                %proposed Ca
    n=X(:,1);                   %proposed spike
    Z=exp(-dC.^2/(2*s^2));      %rho(x)=rho(n)*RHO(C)
    Z(n>0) =Z(n>0)*(1-p);
    Z(n==0)=Z(n==0)*p;
  end
else          %if swarm is defined, use KERNEL DENSITY ESTIMATOR
  if(isempty(X))                  %DRAW SAMPLE FROM rho(X,t)
    %THIS IS THE ONLY PART THAT SHOULD RUN WITH MCMC SAMPLER
    m=size(P.nf{t},1);            %different kernels total    
    ph=P.wf(:,t);                 %weights of different gaussians
    ph=ph/sum(ph);                %normalize, just in case
    phc=cumsum(ph);               %SELECTOR
    
    ind=zeros(N,1);               %draw the kernel from the mixture
    r=rand(N,1);
    for i=1:N ind(i)=find(r(i)<phc,1); end

    n=P.nf{t}(ind,1);             %unwrap the kernel from the mixture 
    C=P.nf{t}(ind,2);
                                  %draw rndized spike-state from kernel,                                  
    n=(rand(N,1)>P.vf(1,t))==(n>0);%TRICK HERE
    C=C+P.vf(2,t)*randn(N,1);     %draw rndized Ca-state from kernel
    Z=[n C];

    %EVALUATE rho(X,t)=P(n)*P(C)
    if(nargout==1) return; end    
    p=sum(ph(P.nf{t}(:,1)==0));  %prob to pull n==0 from diff kernels
                                 %account for switching vf(1,t)
    %P(0)=P(base=0)*P{stay}+P(base=1)*P{switch}
    p=(1-P.vf(1,t))*p+P.vf(1,t)*(1-p);
    
    %prob to pull C from diff Gauss kernels
    s2=P.vf(2,t)^2;    
    A=repmat(C,[1 m]);
    B=repmat(P.nf{t}(:,2),[1 N]);
    A=(A-B').^2/(2*s2);
    A=exp(-A);
    for i=1:m A(:,i)=A(:,i)*ph(i); end%P(kernel)*P(C|kernel)
    p1=sum(A,2);                  %cum prob from all kernels
    
    p1(n>0) =p1(n>0)*(1-p);       %spike state
    p1(n==0)=p1(n==0)*p;
    
    p=p1;
  else                            %CALCULATE rho(X,t)
    n=X(:,1);
    C=X(:,2);
    N=size(C,1);
    
    m=size(P.nf{t},1);            %different kernels total        
    ph=P.wf(:,t);                 %weights of different gaussians
    ph=ph/sum(ph);                %normalize, just in case    
    
    p=sum(ph(P.nf{t}(:,1)==0));  %prob to pull n==0 from diff kernels
                                  %account for switching vf(1,t)
    %P(0)=P(base=0)*P{stay}+P(base=1)*P{switch}
    p=(1-P.vf(1,t))*p+P.vf(1,t)*(1-p);
    
    %prob to pull C from diff Gauss kernels
    s2=P.vf(2,t)^2;    
    A=repmat(C,[1 m]);
    B=repmat(P.nf{t}(:,2),[1 N]);
    A=(A-B').^2/(2*s2);
    A=exp(-A);
    for i=1:m A(:,i)=A(:,i)*ph(i); end%P(kernel)*P(C|kernel)
    p1=sum(A,2);                  %cum prob from all kernels
    
    p1(n>0) =p1(n>0)*(1-p);       %spike state
    p1(n==0)=p1(n==0)*p;
    
    Z=p1;
  end
end


function C = AHill_v1(P,F)
% generalized hill model, inverse
C = (F./(1-F)).^(1/P.n)*P.k_d; C(C<0)  = 0;
