function run_nettime(N,rndinit)
% WRAPPER TIME BIAS script 
if(nargin<1) N=10; end%selection of N
if(nargin<2) rndinit=3711; end%random initialization

fprintf('In run_nettime %i\n',N);
% cd /hmt/sardine/hpc/scratch/stats/users/ym2289/glm-netfit

%%%%%%%%%%  SIMULATION NAME  -- unique to the project
netsim_name_base='data/frame-0717';

%INITIALIZE PARAMETERS
FR=66;              % imaging frame rate
SP=max(0.1,2/N);    % connections sparseness

FRs=[66,1000,5,15,33,100,200];

spkM=1;             % spike-train samples from spk-sampler, for GLM
tmin=1;             % min coupling time-depth, >1
tmax=1;             % max coupling time-depth

Tp=tmax-tmin+1;     % couplings temporal depth

nrange{1}=1:N;      % neurons to process on this proc
id_proc=1;          % this processor
N_proc=1;           % num processors


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SIMULATE DATA
for FR=FRs  
  T=600;              % times to check
  FR_0=FR;
  
  netsim_name=[netsim_name_base,sprintf('-%iFR',FR)];
  fname=[netsim_name,sprintf('-result_%i_%i.mat',N,T)];
  if(exist(fname)==2) continue; end  
  
  fname=[netsim_name,sprintf('-data-%i.mat',N)];        %load sim data, if there
  if(exist(fname)==2) load(fname); else tic;NETFIT_ini;toc; end
  
  T=600;
  FR=FR_0;
  
  
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
  %compute H
  tau_h=0.01;
  for k=1:N
    H{k}(1)=0;
    for t=2:length(H{1})
      H{k}(t)=exp(-1/FR/tau_h)*H{k}(t-1) + n{k}(t);
    end
  end
  
  
  r_full=[];
  weights_iglm=[];            %previous solution carry-over /for ini/
  weights_ispa=[];
  weights_idal=[];
  
  clear lambda dalepr       %RESET FROM PREVIOUS RUN
  lambda=0;
  
  fname=[netsim_name,'temp.mat'];
  ivals=weights_iglm;       %EVALUATE GLM
  bbox=[-10,10];
  NETFIT_GLM_fit
  weights_iglm=ivals;
  Weights_glm=Weights;
  
  result=[];                    % WRITE FULL RESULT
  W=sum(Weights,3);
  c=corrcoef(full(netSim.weights(:)),W(:));
  result.ref=netsim_name;       % reference name
  result.N=N;                   % # neurons
  result.T=netSim.T*netSim.dt;  % train length
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
  
  fname=[netsim_name,sprintf('-result_%i_%i.mat',N,T)];
  save(fname,'result');
end