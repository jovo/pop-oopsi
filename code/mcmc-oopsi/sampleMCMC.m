function [Xr,coll,au]=sampleSMC(P,Y,funcPY,funcPXX,funcG)
%Markov Chain Monte Carlo for conditioned samples in continuous HMM.
% [Xr,Tr,Au]=sampleSMC(p,Yr,@PY,@PXX,@PG) computes chain of states
% conditioned on observations Yr in HMM specified by observation
% probabilities PY and transition probabilities PXX. See help for sample
% spPY and spPXX for definition of the functions setting HMM. Yr is
% cell-array of T observations. "p" is a structure of parameters to be
% passed onto PY and PXX, additionally defining P.iter - number of
% iterations in SMC - and P.grid - number of points to use in stochastic
% "grid". Au is a vector of autocorrelation values, and Tr is a 1xK
% cell-array of K sampled M dimensional states.
% Yuriy Mishchenko 2009 Columbia Un

T=length(Y);
K=P.iter;                     %number of MCMC iterations
N=P.grid;                     %size of stochastic grids

failcnt=0;                    %grid failures counter
failLim=10;                   %grid failures limit

if(nargout>1) coll=cell(1,K); end%sampler history


%GENERATE STOCHASTIC GRID FOR ALL t
X=cell(1,T);                  %form initial stochastic grid /aka pool/
for t=1:T
  [tmp p]=funcG(P,[],t);
  X{t}=zeros(N,length(tmp)+1);%determine state length at t
  X{t}(1,:)=[tmp p];          %add rho(X,t) at the end, for spY
  
  [tmp p]=funcG(P,[],t,N-1);  %generate subsequent X in the grid
  X{t}(2:N,:)=[tmp p];
end


%DRAW HMM SAMPLE ON THE GRID
Xr=[];
it=0;
while(it<K)
  it=it+1;
  Xr_old=Xr;
  
  %draw sample from discretized HMM
  [Xr xr]=sampleMC(P,X,Y,funcPY,funcPXX);
  if(isempty(Xr))             %DEAL WITH GRID PROBLEMS        
    if(failcnt<failLim)       
      fprintf(['Warning: sampleMC returned degenerate sample,',... 
        'will try another grid!\n']);
    else      
      fprintf('Error: could not get nondegenerate sample!\n');
      Xr=[]; coll=[]; au=[]; return
    end
    it=it-1;                  %step back
    failcnt=failcnt+1;        %increment failures
    Xr=Xr_old;                %revert to previous sample    
  else
    failcnt=0;                %fail counter to zero
    if(nargout>1) coll{it}=Xr; end
  end
  
  
  for t=1:T                   %REFORM STOCH GRID /pool/
    k=N;                      %prev chain should be in
    if(~isempty(Xr)) X{t}(1,:)=Xr{t}; k=N-1; end
    [tmp p]=funcG(P,[],t,k);  %form the rest of the grid
    X{t}(N-k+1:end,:)=[tmp p];
  end
  
  if(mod(it,ceil(K/25))==0) fprintf('.'); end
end

%CALCULATE AUTOCOVARIANCE
if(nargout>2) 
  n=length(Xr{1}(:))-1;
  z=zeros(K,n*T);
  for it=1:K
    z1=zeros(T,n);
    for t=1:T z1(t,:)=coll{it}{t}(1:n); end
    z(it,:)=z1(:)';
  end
  z1=xcov(z(:,1),'unbiased'); au=zeros(length(z1),n*T); au(:,1)=z1;
  for it=2:n*T au(:,it)=xcov(z(:,it),'unbiased'); end
end