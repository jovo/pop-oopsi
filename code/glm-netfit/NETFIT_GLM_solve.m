function [Weights,Omega,Rate,J,ivals,lik_tot]=NETFIT_GLM_solve(n,H,range,par,ivals)
%SCRIPT TO SOLVE W GIVEN PARAMETERS
if(nargin<5) ivals=[]; end
    
D=49;                       % coarsening factor, for compression of sample
rate=5;                     % default spontaneous rate
N=length(H);                % # of neurons
T=size(H{1},2);             % available spike history
Tp=par.Tp;                  % coupling temporal depth

%RESET
Rate    = ones(N,1);        % estimated spontaneous rate
Omega   = zeros(N,Tp);      % estimated refractory (self) terms
Weights = zeros(N,N,Tp);    % estimated (cross) coupling terms
J=cell(N,1);                % coupling currents

%%RECOVER WEIGHTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Finding MLE...\n');
warning off all
                            
Sm=1;                       % DESIGN VARIABLES COARSENING/SIZE REDUCTION   
for l=1:N Sm=max(max(H{l}(:)),Sm); end
Sm=Sm/D;                              % bin size 1/D
for l=1:N H{l}=uint8(H{l}/Sm); end    % convert to uint8
if(isfield(par,'lambda')) par.lambda=par.lambda/Sm; end

                            % EXPLANATORY variables
nf=zeros(2+N*par.Tp,T*par.samples,'uint8');
nf(2,:)=uint8(1/Sm);        % spontaneous rate,
for t=par.tmin:par.tmax     % then for all dt-offsets
  for l=1:N
    tmp=[zeros(par.samples,t),H{l}(:,1:end-t)]';
    nf(2+N*(t-par.tmin)+l,:)=tmp(:);  % collapse MxT into 1xM*T
  end  
end

lik_tot=0;
for k=range                           % FOR ALL NEURONS IN RANGE  
  tmp=n{k}';                          %reduce to patterns + frequencies
  nf(1,:)=tmp(:);                     %NETFIT_compress_cols
  [pt,fr]=NETFIT_compress_cols(nf);
  
  bk=zeros(N+1,1); bk(1)=Sm*log(rate);%USE PAST FIT AS INI, if any  
  if(size(ivals,1)>=k) bk=Sm*ivals(k,:)'; end
  
  tic;[bk lik]=NETFIT_glm(pt(1,:),double(pt(2:end,:)),fr,bk,par); toc % MLE
  
  tmp=bk'*double(nf(2:end,:));        % keep track of estimated 
  tmp=reshape(tmp,[T,par.samples]);   % expected value input current
  J{k}=mean(tmp,2)';  

  lik_tot=lik_tot+lik;
  ivals(k,:)=bk;  
  bk=bk'/Sm;                          % scale back the weights
  
  Rate(k) = bk(1);                    % sort stuff into appropriate classes
  for t=par.tmin:par.tmax
    offset=t-par.tmin+1;
    Weights(k,:,offset)=reshape(bk(2+(offset-1)*N:1+offset*N),[1,N,1]);  
    Omega(k,offset)=Weights(k,k,offset);
    Weights(k,k,offset)=0;
  end    
end
warning on all

%clean the mess
clear k l nbuf Hbuf tmp T offset t bk nf nt
    