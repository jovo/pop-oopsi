%RUN MCMC spike train sampler

P=[];
% **************  SAMPLE  *****************  %
T=250;            %spike train length

P.alpha=1;        %fluorescence amplitude
P.beta=0;         %fluorescence background
P.gamma=1e-4;     %fluorescence conversion factor/noise weight
P.zeta=1e-3;      %background noise
P.n=1;            %fluorescence saturation exponent
P.k_d=100;        %fluorescence saturation offset

P.dt=0.05;        %time bin
P.f=30;           %spike generation frequency
P.tau_c=0.2;      %lag of Ca reporter
P.C_0=15;         %Ca background
P.A=30;           %Ca spike height
P.sigma_c=10;     %Ca noise
P.k=log(10);      %log-rate

P.iter=30;        %MCMC iterations
P.grid=25;        %points in stochastic grid
P.minVar=1;       %min kernel dispersion

P.mismatch_penalty=5;%time-error to call a mismatch for bipartite
P.expected_slack=0.5;%fraction of expected slack matches for bipartite

for k=1:1; %trial index
  [Yr,Xr]=prepMC(P,T,@spPY,@spPXX);   %obtain HMM sample

  n=zeros(1,T);                       %unwrap state
  C=zeros(1,T);
  F=zeros(1,T);
  for i=1:T C(i)=Xr{i}(2); n(i)=Xr{i}(1); F(i)=Yr{i}; end;
  figure(k),clf,plot(C/P.A),hold on,stem(n,'r'),plot(F,'g','LineWidth',1.5)

  P.F=F;                              %infer spike train sample
  [Xr Tr Au]=sampleMCMC(P,Yr,@spPY,@spPXX,@spPG);

  n1=zeros(1,length(Xr));             %unwrap inferred state
  C1=zeros(1,length(Xr));
  for t=1:length(Xr) C1(t)=Xr{t}(2); n1(t)=Xr{t}(1); end;
  figure(k),stem(0.9*n1,'c'),plot(C1/P.A,'m:','LineWidth',1.5),pause(0.1);
  legend('true Ca','true spk','true flu','inf spk','inf Ca');

  tgt{k}=n;                           %write down target/estimator
  est{k}=n1;
end

rep=errorsSMC(P,est,tgt);             %determine reconstruction quality
rep.time_inaccuracy=rep.time_inaccuracy*P.dt;

%plot MCMC autocorrelations
figure,plot(sum(Au(:,1:100),2)),title('Neal aucov spikes')
figure,plot(sum(Au(:,101:200),2)),title('Neal aucov calcium')
%display reconstruction valuation
rep