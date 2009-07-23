% this script pulls weights from a set of spike trains
clear, close, clc, fprintf('\nnetSim.WEIGHTS\n')

netsim_name='netSim0315N50S5678mod';
load([netsim_name,'.mat'],'netSim');
load([netsim_name,'INF33Hz.mat'],'Sim','P_infer','n_infer');
FR=Sim.dt; Sim=[];

% %NBAR INFERENCE
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for k=1:length(P_infer) n_infer(k,:)=P_infer{k}.nbar>0.1; end

%DEFINED TEST copy over actual spike trains & check
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% netsim_name='netSim0315N50S5678short';
load([netsim_name,'.mat'],'netSim','n');
T=300;
FR=0.015;
n_infer=zeros(length(n),ceil(T/FR));
for k=1:length(n) 
  t=n{k}*netSim.dt/FR;                  % convert to frame-rate
  t=t(t>0 & t<=size(n_infer,2));        % constrain to frame
  n_infer(k,ceil(t))=1; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Sim.dt          = FR;
Sim.lambda      = 0;
Sim.depth_min   = 0;                    % min coupling time-delay, steps
Sim.depth_max   = 1;                    % max coupling time-delay, steps


%%RECOVER WEIGHTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Finding MLE...\n');
warning off all
Np      = size(n_infer,1);              % number of neurons
mle_spt = zeros(Np,1);                  % spontaneous rate
mle_self= zeros(Np,1);                  % refractory terms
mle_wful= zeros(Np,Np,Sim.depth_max-Sim.depth_min+1);
mle_w   = zeros(Np,Np);                 % weights
flg_start=1;
for k=1:Np                              % for neuron k
  if(exist('tmp.mat') && flg_start)     % load previous data, if any
    flg_start=0;
    load tmp.mat;
    k=k+1;                              % get to next cell
  end
  
  nt=n_infer(k,:);                      % targets
  
  nf=ones(1,size(n_infer,2));           % spontaneous
  for t=Sim.depth_min:Sim.depth_max     % otherwise for offsets
    if(t==0)                            % come up with order in t==0
      order=rand(size(n_infer))<0.5;    % initially random, later Gibbs
      nbuf=n_infer;                     
      nbuf(k,:)=0;                      % no self
      nbuf(order)=0;                    
    else                                % delayed is simple 
      nbuf=[zeros(Np,t),n_infer(:,1:end-t)];
    end
    
    nf=[nf;nbuf];
  end

  tic;bk=NETFIT_glm(nt,nf,[],Sim)';toc  % inference

  mle_spt(k) = bk(1);                   % spontaneous term
  wbuf=zeros(1,Np);                     % coupling terms
  for t=Sim.depth_min:Sim.depth_max
    offset=t-Sim.depth_min+1;
    mle_wful(k,:,offset)=reshape(bk(2+(offset-1)*Np:1+offset*Np),[1,Np,1]);
    wbuf=wbuf+bk(2+(offset-1)*Np:1+offset*Np); 
  end
  mle_self(k)=wbuf(k); wbuf(k)=0;
  mle_w(k,:) = wbuf';                      
  if(mod(k,10)==0) fprintf('='); end
  
  save tmp.mat mle_spt mle_self mle_w mle_wful k
end
warning on all
fprintf('\n');


%%RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure,plot(mle_w(:),netSim.weights(:),'.')
c=corrcoef(mle_w(:),netSim.weights(:)); c=c(3);
title(sprintf('r^2=%.3g',c^2))
fprintf('Simulation result r^2=%.3g\n',c^2);

save([netsim_name,'66HZGLM-base.mat'],'mle_wful','mle_w','mle_spt','mle_self','c');