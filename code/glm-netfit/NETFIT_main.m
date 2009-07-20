%WRAPPER SCRIPT FOR POPULATION MCMC-GIBBS-EM INFERENCE FROM CALCIUM IMAGING
% THIS IS NOT FULLY PARALLELIZED SCRIPT - IT LACKS SYNCHRONIZATION
path(path,'../glm-oopsi/code');%LINK THIS TO PF-filter
path(path,'../mcmc-oopsi');    %LINK THIS to MCMC sampler

%additional solvers, sparse & dale 
if(exist('rate_b')~=1) rate_b=5; end        %expected spike rate
if(exist('flgSparse')~=1) flgSparse=0; end  %find sparse solution
if(exist('flgDale')~=1) flgDale=0; end      %find Dale solution
if(exist('flgData')~=1) flgData=0; end      %default is simulation

if(flgData==0)              % use to simulate dataset
  fname=[netsim_name,'-data.mat'];
  %load sim data, if there, or simulate dataset -- one node only
  if(id_proc==1 && exist(fname)~=2) tic; NETFIT_ini; toc; end
  while(exist(fname)~=2) fprintf('waiting...\n'); pause(15); end
  load(fname);
end


%INNER VARIABLES
n=cell(N,1);                % spikes array
H=cell(N,1);                % spike history array
J=cell(N,1);                % coupling currents
M=cell(N,1);                % estimated calcium models
Rate    = ones(N,1);        % estimated spontaneous rate
Omega   = zeros(N,Tp);      % estimated refractory (self) terms
Weights = zeros(N,N,Tp);    % estimated (cross) coupling terms

Weights_glm = zeros(N,N,Tp);% estimated (cross) coupling terms
Weights_spa = zeros(N,N,Tp);% estimated (cross) coupling terms, base glm
Weights_dal = zeros(N,N,Tp);% estimated (cross) coupling terms, base glm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cnt=0;                      % loops counter
iflg=1;                     % main-EM-loop flag
weights_iglm=[];            % previous solution carry-over /for solver ini/
weights_ispa=[];
weights_idal=[];
Weights_old=repmat(Inf,size(Weights));
while(iflg)         
  if(~strcmp(mode,'base'))
                            % PERFORM ESATIMATION OF M-MODELS
                            % incoming fluorescence traces in C{k}
                            % incoming currents in J{k}
                            % outgoing estimated calcium models M{k}
    fname=[netsim_name,sprintf('-models%iproc%i.mat',cnt,id_proc)];
    if(exist(fname)==2) load(fname); else NETFIT_EM_loop_M; end
  
    syncdata={'M'};           % SYNCHRONIZE
    NETFIT_sync
  end

                            % DRAW SAMPLE SPIKE TRAINS
                            % incoming fluorescence traces in C{k}                     
                            % incoming currents in J{k}
                            % outgoing sample of spike trains n{k} ----|
                            % outgoing sample of spike histories H{k}--|
  fname=[netsim_name,sprintf('-sample_%s%iproc%i.mat',mode,cnt,id_proc)];
  if(exist(fname)==2) load(fname); else NETFIT_spk_sampler; end
  
  syncdata={'H','n'};       % SYNCHRONIZE
  NETFIT_sync
  
  
  clear lambda dalepr       %RESET FROM PREVIOUS RUN  
  
  
                            % EVALUATE GLM
                            % ALSO evaluate coupling currents J
                            % incoming sample of spike trains n{k} ----|
                            % incoming sample of spike histories H{k}--|
                            % outgoing couping weights W, rates R, refr Omega
                            % outgoing coupling currents J  
  ivals=weights_iglm;       
  bbox=[-10,15];            % bounding box on weigths
  fname=[netsim_name,sprintf('-glm_%s%iproc%i.mat',mode,cnt,id_proc)];
  if(exist(fname)==2) load(fname); else NETFIT_GLM_fit; end  
  
  syncdata={'W','V'};       % SYNCHRONIZE
  W=cell(N,1); V=cell(N,1); % which variables will synch
  for k=nrange{id_proc} W{k}=Weights(k,:); V{k}=ivals(k,:); end
  NETFIT_sync               % sync here
  for k=1:N Weights(k,:)=W{k}; ivals(k,:)=V{k}; end  
  weights_iglm=ivals;       % unwrap synced variables
  Weights_glm=Weights;
    
            
  if(flgSparse)             % EVALUATE GLM WITH SPARSE PRIOR
    if(exist('lambda_b')==1) lambda=lambda_b; else lambda=20; end

    % NETFIT_GLM_lambda       % find/use best lambda
    % ivals=weights_ispa;     % if commented, use previous solution as init
    bbox=[];               % bounding box for weights
    fname=[netsim_name,sprintf('-spa_%s%iproc%i.mat',mode,cnt,id_proc)];
    if(exist(fname)==2) load(fname); else NETFIT_GLM_fit; end
    
    syncdata={'W','V'};     % SYNCHRONIZE
    W=cell(N,1); V=cell(N,1);
    for k=nrange{id_proc} W{k}=Weights(k,:); V{k}=ivals(k,:); end
    NETFIT_sync
    for k=1:N Weights(k,:)=W{k}; ivals(k,:)=V{k}; end
    weights_ispa=ivals;
    Weights_spa=Weights;
  else
    Weights_spa=Weights;
  end
  
  if(flgDale)
    NETFIT_GLM_signs       % EVALUATE GLM WITH DALE PRIOR
    % ivals=weights_idal;    % if commented, use previous solution as init
    bbox=[];               % bounding box for weights    
    fname=[netsim_name,sprintf('-dal_%s%iproc%i.mat',mode,cnt,id_proc)];
    if(exist(fname)==2) load(fname); else NETFIT_GLM_fit; end
    
    syncdata={'W','V','J'}; % SYNCHRONIZE
    W=cell(N,1); V=cell(N,1);
    for k=nrange{id_proc} W{k}=Weights(k,:); V{k}=ivals(k,:); end
    NETFIT_sync
    for k=1:N Weights(k,:)=W{k}; ivals(k,:)=V{k}; end
    weights_idal=ivals;
    Weights_dal=Weights;
  else
    Weights_dal=Weights;
  end

  NETFIT_eval_term          % evaluate termination conditions
                            % outgoing outer EM-loop flag flg                            
end

if(id_proc>1) return; end  % ONLY ONE NODE NEED TO PREP THAT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
result=[];                  % WRITE FULL RESULT
result.ref=netsim_name;     % reference name
result.F=F;                 % fluorescence data
result.GT=n_GT;             % spikes ground truth
result.GT_K=netSim.K;       % downsampling factor in n_GT
result.GT_W=netSim.weights; % true weights

result.M=M;                 % infered calcium models  
result.n=n;                 % infered spike trains
result.H=H;                 % infered spike histories

result.Weights_glm=Weights_glm; % infered connectivity weights
result.Weights_spa=Weights_spa; % infered connectivity weights
result.Weights_dal=Weights_dal; % infered connectivity weights

W=sum(result.Weights_dal,3);
c=corrcoef(full(netSim.weights(:)),W(:));
result.Q=c(3)^2;            % corr-coef-square

vname=['result_',mode];
eval([vname,'=result;']);

fname=[netsim_name,sprintf('-result_%s.mat',mode)];
save(fname,vname);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
result=[];                  % WRITE COMPACT RESULT
result.ref=netsim_name;     % reference name
result.GT_W=netSim.weights; % true weights
result.Weights_glm=Weights_glm; % infered connectivity weights
result.Weights_spa=Weights_spa; % infered connectivity weights
result.Weights_dal=Weights_dal; % infered connectivity weights
result.Q=c(3)^2;            % corr-coef-square
 
vname=['result_',mode];
eval([vname,'=result;']);

fname=[netsim_name,'-result_bf.mat'];
if(exist(fname)~=2) save(fname,vname); else save(fname,'-append',vname); end