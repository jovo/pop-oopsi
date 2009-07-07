function [bk lik_r]=NETFIT_glm(n,H,fr,x0,par)
%Function fits GLM to spiking data.
% [bk lik_r]=NETFIT_glm(n,H[,fr,x0,options]) fits GLM with parameters bk 
% to explain spike train n with spike histories H. n is 1xT vector and 
% H is NxT matrix; x0 may be passed as initial guess for bk; fr is 1xT
% vector of frequencies of stat patterns in [n;H].
% options is parameters structure, containing following fields: 
% dt - the time step, lambda - L_1 sparsifying weight.
%
% Note, in order to look at M samples, stack samples one after another 
% into 1xT*M vector n and NxT*M matrix H.
%
% Y. Mishchenko Columbia University Department of Statistics 2009

N=size(H,1);                          % # of explanatory variables
T=size(H,2);                          % # of samples

if(nargin<5) par=[]; end              % initialize params
if(~isfield(par,'dt')) par.dt=0.01; end
if(~isfield(par,'lambda')) par.lambda=0; end% L1 sparse prior
if(~isfield(par,'signs')) par.signs=[]; end % dale prior, row signs
if(~isfield(par,'samples')) par.samples=1; end %number of samples
if(nargin<3 || isempty(fr)) fr=ones(1,T); end %default frequencies
if(nargin>=4 && ~isempty(x0)) w=x0; else w=rand(N,1); end% initialize weights

if(par.lambda==0)
  options = optimset('Display','off','GradObj','off','Hessian','off',...
   'TolFun',1e-6,'MaxIter',10000,'MaxFunEvals',10000); % pass some options to solver
else
  options = optimset('Display','off','GradObj','on','Hessian','off',...
   'TolFun',1e-6,'MaxIter',10000,'MaxFunEvals',10000);% pass some options to solver
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SPARSE/DALE (CON) & NONSPARSE (UNC) PROBLEMS ARE SOLVED SEPARATELY
% nonsparse problem is "min -L(w)"
% sparse (L1) problem is "min -L(w)+Lambda*|w|", we approach it as
% more regular "min -L(w)+Lambda*sum t_i" s.t. "w_i-t_i<=0, -w_i-t_i<=0"
cnt=0;
if(par.lambda==0 && isempty(par.signs~=0))% VANILLA MLE
  %[bk lik_r]  = fminunc(@fbk,w,options); %this crashed in strong coupling regime
  [bk lik_r]=fmincon(@fbk,w,[],[],[],[],-repmat(15,size(w)),repmat(15,size(w)),[],options);
elseif(par.lambda==0)                 % DALE PRIOR ALONE
  d=zeros(1,N-1);                     % initialize constrained-min problem
  d(find(par.signs<0))=1;             % dale's w<0 sub  w<0
  d(find(par.signs>0))=-1;            % dale's w>0 sub -w<0
  d=[-1,d];                           % rate>0
  A=sparse(diag(d));          
  a=max(abs(A),[],2);
  A=A(a>0,:);                         % only keep constrained vars
  
  b=zeros(size(A,1),1);               % constraints RHS
  
                                      % feasible initial w
  w(find(par.signs>0)+1)=max(0,w(find(par.signs>0)+1));
  w(find(par.signs<0)+1)=min(0,w(find(par.signs<0)+1));
  
  [bk lik_r]  = fmincon(@fbk,w,A,b,[],[],[],[],[],options);  
else                                  % L1 PRIOR + DALE PRIOR
  d=sparse(diag(ones(1,N-1)));        % initialize constrained-min problem
  d0=d; d0(:,find(par.signs<0))=0;    % dale's w<0 sub  w-t<0 with  w<0
  d1=d; d1(:,find(par.signs>0))=0;    % dale's w>0 sub -w-t<0 with -w<0
  d=[zeros(N-1,1),d];                 % rate is not regularized, w(1)
  
  A=[[d,-d0];[-d,-d1]];               % matrix of constraints
  A=[zeros(1,2*N-1);A];               % constrain the rate to be positive
  A(1,1)=-1;                          % -w(1)<0 i.e. -rate<0
  
  b=zeros(2*N-1,1);                   % constraints RHS
  
                                      % feasible initial [w,t]   
  w(find(par.signs>0)+1)=max(0,w(find(par.signs>0)+1));
  w(find(par.signs<0)+1)=min(0,w(find(par.signs<0)+1));
  w=[w;abs(w(2:end))+0.1/par.lambda]; 
  
  [bk lik_r]  = fmincon(@fbk,w,A,b,[],[],[],[],[],options);
  
  bk=bk(1:N);                         % drop dummy t's
end
fprintf('cost %i',cnt);
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CALCULATING LIKELIHOOD
function [lik dlik] = fbk(w)    
  cnt=cnt+1;
  eps=1e-6;                           % safety margin
  
  %parse incoming point; if L1 - spit into w & t
  if(par.lambda>0) w0=w(1:N); t0=w(N+1:end); else w0=w; end
  
  s=w0'*H;                            % LM
  log_pr_nosp=-exp(s)*par.dt;         % log-Pr no spike
                                      % WE REALLY DON'T CARE
                                      % ABOUT FOLLOWING AT n==0
  pr_nosp=exp(log_pr_nosp(n>0));      % Pr no spike
  pr_sp=max(eps,1-pr_nosp);           % Pr spike
  log_pr_sp=log(pr_sp);               % log Pr spike

  lik=sum(log_pr_sp.*fr(n>0))+sum(log_pr_nosp(n==0).*fr(n==0));%lik      
  lik=-lik/par.samples;               % minimization, adj samples # (!)
  
  if(par.lambda>0) lik=lik+par.lambda*sum(t0); end%L1-term (!)
  
  if(nargout>1)                       % gradient
    term_sp=repmat(1./pr_sp.*pr_nosp.*(-log_pr_nosp(n>0)).*fr(n>0),[N 1]).*H(:,n>0);
    term_sp=sum(term_sp,2);           % diff of lik-spikes
    
    term_nosp=repmat(log_pr_nosp(n==0).*fr(n==0),[N 1]).*H(:,n==0);
    term_nosp=sum(term_nosp,2);       % diff of lik-nospikes
    dlik=-(term_sp+term_nosp)/par.samples;
    
    if(par.lambda>0) dlik=[dlik;repmat(par.lambda,size(t0))]; end%t-terms
  end
end
    
end
