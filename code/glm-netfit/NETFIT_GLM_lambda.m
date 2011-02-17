%SCRIPT TO CHOOSE BEST LAMBDA FOR GIVEN W USING GT

par=[];
par.lambda=0;                     % L1-regularizer strength
par.tmin=tmin;                    % min coupling time-depth
par.tmax=tmax;                    % max coupling time-depth
par.dt=netSim.FR;                 % time-step
par.samples=spkM;                 % samples in the data
par.Tp=Tp;                        % coupling temporal depth

Ntest=10;                         % reduced # of rows to go through
lambdas=[0,5,10,20,30,40,50,75,100];% L1-weights to check
lambdas=[lambdas,1/mean(netSim.weights(:)~=0)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Finding Best Lambda...\n');
cc=zeros(size(lambdas));          % obtained corr. coefs.
for l=1:length(lambdas);
  par.lambda=lambdas(l);          % PICK WEIGHT
  
  if(par.lambda==0 && exist('Weights_glm')==1)
    W=Weights_glm;
  else                            % solve only for Ntest neurons ;)
    W=NETFIT_GLM_solve(n,H,1:Ntest,par,[]);
  end 

  W=W(1:Ntest,:);  
  tmp=netSim.weights(1:Ntest,:);  % get the fragment of true matrix
  c=corrcoef(W(:),tmp(:));        % get corr. coef. for this lambda
  cc(l)=c(3);
  fprintf('C.C. %g\n',c(3));
end

[tmp,k]=max(cc);                    % pick best weight
lambda=lambdas(k);
fprintf('Best lambda %g\n',lambda);

%clean the mess
clear Ntest lambdas cc tmp W lc par
