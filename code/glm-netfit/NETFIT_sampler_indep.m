%SCRIPT TO DRAW SAMPLE OF SPIKE HISTORIES FOR ESTIMATING W
% --- factorized/independent approximation, using PF

Sim=[];
Sim.N       = 50;                       % # of particles
Sim.dt      = netSim.FR;                % time step size
Sim.freq    = 1;                        % # of time steps between observations
Sim.M       = 0;                        % number of spike history dimensions
Sim.pf      = 1;                        % use conditional sampler (not prior) when possible

                                        % don't infer anything here
Sim.Mstep   = 0;                        % do M-step
Sim.C_params= 0;                        % whether to estimate calcium parameters {tau,A,C_0,sig}
Sim.n_params= 0;                        % whether to estimate rate governing parameters {b,k}
Sim.h_params= 0;                        % whether to estimate spike history parameters {h}
Sim.F_params= 0;                        % whether to estimate observation parameters {alpha,beta}
Sim.G_params= 0;                        % whether to estimate observation parameters {gamma}
Sim.MaxIter = 1;                        % max # of EM iterartions

Sim.Scan    = 1;                        % end-of-frame data

for k=nrange{id_proc}
  %% 2) initialize parameters
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  P=M{k};
  if(isfield(P,'omega')) Sim.M=1; end
  if(isfield(P,'FQ')) Sim.freq=P.FQ; Sim.dt=netSim.FR/P.FQ; end
  
  
  Ft=repmat(F{k}(:)',[Sim.freq 1]); Ft=Ft(:);
  Sim.T       = length(Ft);               % # of samples
  Sim.Nsec    = Sim.T*Sim.dt;             % # of actual seconds
  Sim.T_o     = Sim.T;                    % # of observations
  Sim.tvec    = Sim.dt:Sim.dt:Sim.Nsec;   % vector of times

  Sim.n=repmat(NaN,[1,Sim.T]);            % for plotting purposes in PF
  if(exist('n_GT')==1) Sim.n(min(Sim.T,ceil(n_GT{k}/netSim.K*Sim.freq)))=1; end

  Sim.StimDim=1;  
  Sim.x=ones(1,Sim.T);
  if(~isempty(J{k})) Sim.x=J{k}(1:Sim.T); end % past-estimated incoming currents
    

  %% 3) marginal sample using PF
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  warning off
  [I.M I.P]   = GOOPSI_main_v1_0(Ft(1:Sim.T),P,Sim);
  fprintf('\n');
  warning on
  
  %% 4) choose sample for joints
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % factorized/independent approximation means mix & match
  % choose samples
  tab=cumsum(I.M.w,1);                    % cum sums to choose from
  H_loc=zeros(spkM,Sim.T);
  n_loc=false(spkM,Sim.T);
  for l=1:spkM
    r=rand(1,Sim.T);
    tmp=repmat(1:Sim.N,[Sim.T 1])';       % indices template
    tmp(repmat(r,[Sim.N 1])>tab)=Inf;     % Inf indices above rand
    ind=min(tmp,[],1);                    % pick smallest remaining index
    ind=sub2ind(size(I.M.w),ind,1:Sim.T); % convert to 2D indices

    H_loc(l,:)=I.M.h(ind);                % our sample, done
    n_loc(l,:)=I.M.n(ind);
  end
  n{k}=n_loc;  
  H{k}=[H_loc(:,2:end),zeros(spkM,1)];    % sample of histories for k - MxT array    
                                          % shift forward by 1 since 
                                          % NETFIT_GLM_fit will shift by 1
end

save(fname,'n','H');

%clean the mess
clear P Sim I tab H_loc n_loc k l r tmp ind
