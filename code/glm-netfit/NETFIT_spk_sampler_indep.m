%SCRIPT TO DRAW SAMPLE OF SPIKE HISTORIES FOR ESTIMATING W
% --- factorized/independent approximation, using PF

fname=[netsim_name,sprintf('-SMPL-%i.mat',cnt)];
if(exist(fname)==2) load(fname); return; end

Sim=[];
Sim.dt      = netSim.FR;                % time step size
Sim.freq    = 1;                        % # of time steps between observations
Sim.N       = 50;                       % # of particles
Sim.M       = 1;                        % number of spike history dimensions
Sim.pf      = 1;                        % use conditional sampler (not prior) when possible

                                        % don't infer anything here
Sim.Mstep   = 0;                        % do M-step
Sim.C_params= 0;                        % whether to estimate calcium parameters {tau,A,C_0,sig}
Sim.n_params= 0;                        % whether to estimate rate governing parameters {b,k}
Sim.h_params= 0;                        % whether to estimate spike history parameters {h}
Sim.F_params= 0;                        % whether to estimate observation parameters {alpha,beta}
Sim.G_params= 0;                        % whether to estimate observation parameters {gamma}
Sim.MaxIter = 1;                       % max # of EM iterartions

Sim.Scan    = 1;                        % end-of-frame data

for k=1:N
  %ESTIMATE CELL MODEL FOR NEURON k
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  Sim.T       = length(F{k});             % # of samples
  Sim.Nsec    = Sim.T*Sim.dt;             % # of actual seconds
  Sim.T_o     = Sim.T;                    % # of observations
  Sim.tvec    = Sim.dt:Sim.dt:Sim.Nsec;   % vector of times

  Sim.n=repmat(NaN,[1,Sim.T]);            % for plotting purposes in PF
  if(exist('n_GT')==1) Sim.n(min(Sim.T,ceil(n_GT{k}/netSim.K)))=1; end

  if(isempty(J{k}))                       % past-estimated incoming currents
    Sim.x=ones(1,Sim.T);
    Sim.StimDim=1;
  else
    Sim.x=[ones(1,Sim.T);J{k}];
    Sim.StimDim=size(Sim.x,1);
  end
  
  %% 2) initialize parameters
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  P=M{k};

  %% 3) estimate parameters
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  warning off
  [I.M I.P]   = GOOPSI_main_v1_0(F{k},P,Sim);
  fprintf('\n');
  warning on
  
  %% 4) choose sample
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % factorized/independent approximation means mix & match
  % choose samples
  tab=cumsum(I.M.w,1);                    % cum sums to choose from
  out_H=zeros(spkM,Sim.T);
  out_n=false(spkM,Sim.T);
  for l=1:spkM
    r=rand(1,Sim.T);
    tmp=repmat(1:Sim.N,[Sim.T 1])';       % indices template
    tmp(repmat(r,[Sim.N 1])>tab)=Inf;     % Inf indices above rand
    ind=min(tmp,[],1);                    % pick smallest remaining index
    ind=sub2ind(size(I.M.w),ind,1:Sim.T); % convert to 2D indices

    out_H(l,:)=I.M.h(ind);                % our sample, done
    out_n(l,:)=I.M.n(ind);
  end
  n{k}=out_n;  
  H{k}=[out_H(:,2:end),zeros(spkM,1)];    % sample of histories for k - MxT array    
                                          % shift forward by 1 since 
                                          % NETFIT_GLM_fit will shift by 1
end

save(fname,'n','H');

%clean the mess
clear P Sim I tab out_H out_n k l r tmp ind
