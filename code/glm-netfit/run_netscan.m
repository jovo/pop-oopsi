function run_netscan(N,rndinit)
%LAUNCH script for netfit package-wrapper NETFIT_main
%THIS IS NOT MULTIPROC SCRIPT, SEE run_netprobe FOR MULTIPROC
if(nargin<1) N=10; end%selection of N
if(nargin<2) rndinit=3711; end%random initialization

fprintf('In run_netprobe %i\n',N);
% cd /hmt/sardine/hpc/scratch/stats/users/ym2289/glm-netfit

%%%%%%%%%%  SIMULATION NAME  -- unique to the project
netsim_name='data/scan-0622';

%INITIALIZE PARAMETERS
FR=66;              % imaging frame rate
FR=33;        
SP=max(0.1,2/N);    % connections sparseness

T_range=[300:300:3600]; % times to check

spkM=1;             % spike-train samples from spk-sampler, for GLM
tmin=1;             % min coupling time-depth, >1
tmax=1;             % max coupling time-depth

Tp=tmax-tmin+1;     % couplings temporal depth

nrange{1}=1:N;      % neurons to process on this proc
id_proc=1;          % this processor
N_proc=1;           % num processors


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SIMULATE DATA
T=max(T_range);
fname=[netsim_name,sprintf('-data-%i.mat',N)];        %load sim data, if there
if(exist(fname)==2) load(fname); else tic;NETFIT_ini;toc; end


%INNER VARIABLES
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
%DOWNSAMPLE SPIKES
cnt=0;                      % for compatibility purpose                            
fname=[netsim_name,sprintf('-sample-%i.mat',N)];
if(exist(fname)==2) load(fname); else NETFIT_sampler_base; end 

H_full=H;
n_full=n;
r_full=[];
weights_iglm=[];            %previous solution carry-over /for ini/
weights_ispa=[];
weights_idal=[];
for Ti=T_range
  for k=1:N                 %DEFINE SPK SAMPLE FOR Ti
    H{k}=H_full{k}(1:round(Ti/netSim.FR));
    n{k}=n_full{k}(1:round(Ti/netSim.FR));
  end
  
  clear lambda dalepr       %RESET FROM PREVIOUS RUN
  
  fname=[netsim_name,'temp.mat'];
  ivals=weights_iglm;       %EVALUATE GLM 
  NETFIT_GLM_fit      
  weights_iglm=ivals;  
  Weights_glm=Weights;

  NETFIT_GLM_lambda         %EVALUATE GLM-SPARSE, BEST lambda
  ivals=weights_ispa;
  NETFIT_GLM_fit
  weights_ispa=ivals;
  Weights_spa=Weights;  
  
  NETFIT_GLM_signs          %EVALUATE GLM-DALE  
  ivals=weights_idal;
  NETFIT_GLM_fit
  weights_idal=ivals;
  Weights_dal=Weights;  
  
  result=[];                    % WRITE FULL RESULT
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

vname=sprintf('result_%i',N); % UPDATE FILE
eval([vname,'=r_full;']);
fname=[netsim_name,'-probe.mat'];
if(exist(fname)~=2) save(fname,vname); else save(fname,'-append',vname); end
