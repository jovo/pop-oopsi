function [Xr,coll,au]=sampleSMC(P,Y,funcPY,funcPXX,funcG)
%Sequential Monte Carlo for conditioned samples in continuous HMM.
%[Xr,Tr,Au]=sampleSMC(p,Yr,@PY,@PXX,@PG) computes chain of states 
%conditioned on observations Yr in HMM specified by observation 
%probabilities PY and transition probabilities PXX. See help for sample 
%spPY and spPXX for definition of the functions setting HMM. Yr is 
%cell-array of T observations. "p" is a structure of parameters to be 
%passed onto PY and PXX, additionally defining P.iter - number of 
%iterations in SMC - and P.grid - number of points to use in stochastic 
%"grid". Au is a vector of autocorrelation values, and Tr is a 1xK
%cell-array of K sampled M dimensional states. 
%Yuriy Mishchenko 2009 Columbia Un

T=length(Y);
K=P.iter;                     %number of iterations of MC till convergence
N=P.grid;                     %size of stochastic grid to use

rejcnt=0;                     %grid rejections counter                        
failcnt=0;                    %grid failures counter
failLim=10;                   %grid failures limit
if(nargout>1) coll=cell(1,K); end%sampler history, as needed

X=cell(1,T);                  %form initial stochastic grid /aka pool/
for t=1:T 
  [tmp p]=funcG(P,[],t);
  X{t}=zeros(N,length(tmp)+1);%determine state length 
  X{t}(1,:)=[tmp p];          %add rho(X,t) at the end, for spY
  
  for k=2:N                   %generate subsequent X in the grid
    [tmp p]=funcG(P,[],t); 
    X{t}(k,:)=[tmp p];  
  end
end

Xr=[];
for it=1:K
  Xr_old=Xr;
  [Xr xr]=sampleMC(P,X,Y,funcPY,funcPXX);%draw sample MC
  if(isempty(Xr))             %if here, means problems with sampleMC
    it=it-1;                  %step back in "it
    failcnt=failcnt+1;        %increment failures counter    
    Xr=Xr_old;                %revert to previous sample chain
    
    if(failcnt<failLim)
      fprintf('Warning: sampleMC returned degenerate sample, will try another grid!\n'); 
    else
      Xr=[]; coll=[]; au=[];      
      fprintf('Error: could not get nondegenerate sample!\n');
      return
    end
  else 
    failcnt=0;                %fail counter to zero
    if(it>1 && max(xr)==1) rejcnt=rejcnt+1; end%record rejections   
    if(nargout>1) coll{it}=Xr; end
  end
  
  if(it<K)
    for t=1:T                 %reform stochastic grid /pool/
      if(~isempty(Xr))        
        ko=2;                 
        X{t}(1,:)=Xr{t};      %prev chain should be in
      else ko=1; end
      for k=ko:N              %form rest of the grid
        [tmp p]=funcG(P,[],t);
        X{t}(k,:)=[tmp p];
      end
    end
    if(mod(it,ceil(K/25))==0) fprintf('.'); end
  else
    fprintf('a.r.=%.3g\n',1-rejcnt/(K-1));%print rejection rate
  end
end

if(nargout>2) %obtain autocovariance
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