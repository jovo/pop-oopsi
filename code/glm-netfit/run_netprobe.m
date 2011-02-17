function run_netprobe(N,id_proc,N_proc,rndinit)
% WRAPPER LAUNCH script for netfit package-wrapper NETFIT_main
if(nargin<1) N=10; end%selection of N
if(nargin<2) id_proc=1; end %processor to execute
if(nargin<3) N_proc=1; end  % # of processors to use
if(nargin<4) rndinit=3711; end%random initialization

fprintf('In run_netprobe %i\n',N);
% cd /hmt/sardine/hpc/scratch/stats/users/ym2289/glm-netfit

%%%%%%%%%%  SIMULATION NAME  -- unique to the project
netsim_name='data/probe-0622';
sync_name=[netsim_name,sprintf('-%i-sync',N)];

%INITIALIZE PARAMETERS
FR=66;              % imaging frame rate
SP=max(0.1,2/N);    % connections sparseness

T_range=[300,600,1800,3600]; % times to check

spkM=1;             % spike-train samples from spk-sampler, for GLM
tmin=1;             % min coupling time-depth, >1
tmax=1;             % max coupling time-depth

Tp=tmax-tmin+1;     % couplings temporal depth

D=ceil(N/N_proc);
for i=1:N_proc nrange{i}=(1+(i-1)*D):min(N,i*D); end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIMULATE DATA -- one node
T=max(T_range);     
fname=[netsim_name,sprintf('-data-%i.mat',N)];  %sim data
if(id_proc==1 && exist(fname)~=2) tic; NETFIT_ini; toc; end
while(exist(fname)~=2) pause(15); end
load(fname);

                            % INNER VARIABLES
n=cell(N,1);                % spikes array
H=cell(N,1);                % spike history array
M=cell(N,1);
J=cell(N,1);                % coupling currents
Rate    = ones(N,1);        % estimated spontaneous rate
Omega   = zeros(N,Tp);      % estimated refractory (self) terms
Weights = zeros(N,N,Tp);    % estimated (cross) coupling terms

Weights_glm = zeros(N,N,Tp);% estimated (cross) coupling terms
Weights_spa = zeros(N,N,Tp);% estimated (cross) coupling terms, base glm
Weights_dal = zeros(N,N,Tp);% estimated (cross) coupling terms, base glm


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DOWNSAMPLE SPIKES -- one node
cnt=0;                      % for compatibility purpose                            
fname=[netsim_name,sprintf('-sample-%i.mat',N)];
if(id_proc==1 && exist(fname)~=2) 
  tmp=nrange{1};
  nrange{1}=1:N;
  NETFIT_sampler_base;
  nrange{1}=tmp;
end
while(exist(fname)~=2) pause(15); end
load(fname);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H_full=H;
n_full=n;
r_full=[];
weights_iglm=[];            % previous solution carry-over /for ini/
weights_ispa=[];
weights_idal=[];
for Ti=T_range
  for k=1:N                 % FORM SPK SAMPLE FOR Ti
    H{k}=H_full{k}(1:round(Ti/netSim.FR));
    n{k}=n_full{k}(1:round(Ti/netSim.FR));
  end
  
  
  clear lambda dalepr       % RESET FROM PREVIOUS RUN
  
  
  fname=[netsim_name,sprintf('-glm_%i_%iproc%i.mat',N,Ti,id_proc)];
  ivals=weights_iglm;       
  if(exist(fname)==2) load(fname); else NETFIT_GLM_fit; end  
  
  syncdata={'W','V'};       % SYNCHRONIZE
  W=cell(N,1); V=cell(N,1);
  for k=nrange{id_proc} W{k}=Weights(k,:); V{k}=ivals(k,:); end
  NETFIT_sync
  for k=1:N Weights(k,:)=W{k}; ivals(k,:)=V{k}; end  
  weights_iglm=ivals;       
  Weights_glm=Weights;

  
  NETFIT_GLM_lambda         % EVALUATE GLM-SPARSE, SET lambda
  fname=[netsim_name,sprintf('-spa_%i_%iproc%i.mat',N,Ti,id_proc)];
  ivals=weights_ispa;       
  if(exist(fname)==2) load(fname); else NETFIT_GLM_fit; end  
  
  syncdata={'W','V'};       % SYNCHRONIZE
  W=cell(N,1); V=cell(N,1);
  for k=nrange{id_proc} W{k}=Weights(k,:); V{k}=ivals(k,:); end
  NETFIT_sync
  for k=1:N Weights(k,:)=W{k}; ivals(k,:)=V{k}; end    
  weights_ispa=ivals;
  Weights_spa=Weights;  
  
  
  NETFIT_GLM_signs          % EVALUATE GLM-DALE  
  fname=[netsim_name,sprintf('-dal_%i_%iproc%i.mat',N,Ti,id_proc)];
  ivals=weights_idal;       
  if(exist(fname)==2) load(fname); else NETFIT_GLM_fit; end  
  
  syncdata={'W','V'};       % SYNCHRONIZE
  W=cell(N,1); V=cell(N,1);
  for k=nrange{id_proc} W{k}=Weights(k,:); V{k}=ivals(k,:); end
  NETFIT_sync
  for k=1:N Weights(k,:)=W{k}; ivals(k,:)=V{k}; end    
  weights_idal=ivals;  
  Weights_dal=Weights;  
  
  if(id_proc>1) continue; end

  result=[];                    % WRITE FULL RESULT, ONE NODE
  W=sum(Weights,3);
  c=corrcoef(full(netSim.weights(:)),W(:));
  result.ref=netsim_name;       % reference name
  result.N=N;                   % # neurons
  result.T=Ti;                  % train length
  result.lambda=lambda;         % best L1-weight
  result.Q=c(3)^2;              % fin corr-coef-square
  result.GT=n_GT;               % spikes ground truth
  result.GT_K=netSim.K;         % downsampling factor in n_GT
  result.GT_W=netSim.weights;   % true weights  
  result.n=n_full;              % downsampled spike trains
  result.H=H_full;              % downsampled spike trains
  result.Weights_glm=Weights_glm; % infered connectivity weights  
  result.Weights_spa=Weights_spa; % infered connectivity weights  
  result.Weights_dal=Weights_dal; % infered connectivity weights    
  
  vname=sprintf('result_%i',Ti);  % UPDATE FILE
  eval([vname,'=result;']);  
  fname=[netsim_name,sprintf('-result_%i_%i.mat',N,Ti)];
  save(fname,vname);
  

  result=[];                    % WRITE COMPACT RESULT
  result.ref=netsim_name;       % reference name
  result.N=N;                   % # neurons
  result.T=Ti;                  % train length
  result.lambda=lambda;         % best L1-weight  
  result.Q=c(3)^2;              % corr-coef-square
  result.GT_W=netSim.weights;   % true weights  
  result.Weights_glm=Weights_glm; % infered connectivity weights  
  result.Weights_spa=Weights_spa; % infered connectivity weights  
  result.Weights_dal=Weights_dal; % infered connectivity weights    
  
  if(isempty(r_full)) r_full=result; else r_full(end+1)=result; end
end

if(id_proc>1) return; end

vname=sprintf('result_%i',N);   % UPDATE FILE
eval([vname,'=r_full;']);
fname=[netsim_name,'-probe.mat'];
if(exist(fname)~=2) save(fname,vname); else save(fname,'-append',vname); end
