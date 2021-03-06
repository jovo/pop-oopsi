%% basic setups
aim=0;              %set this to 1 if aiming cell params

N=50;               %# neurons
sp=0.1;             %sparseness
T=60;               %sim length, sec
flags=false(1,N);   %which are excitatory
flags(1:ceil(0.8*N))=1;

dt=0.1;             %Euler time step, ms

if(aim) T=1; sp=0; end %if aiming...
%% define connectivity
W=exprnd(1,N,N);    %exp weights
W(rand(N,N)>sp)=0;  %make sparse
for k=1:N W(k,k)=0; end

%% IF parameters
tauS=1;         %time constant for synapse alpha function, ms
VE=40;          %exc reversal potential, mV
scaleE=0.315;    %excitation strength [a.u.]

VI=-90;         %inh reversal potential, mV
scaleI=1;       %inhibition strength [a.u.]

tauR=10;        %time constant membrane relaxation, ms
VR=-65;         %membrane rest potential

Vthr=-35;       %spiking threshold
Vres=-65;         %spike reset potential, mV


poisE=1000;     %excitatory Poisson input, Hz
poisI=300;      %inhibitory Poisson input, Hz

%% some simple checks & auxiliary stuff
dttrue=dt/1000;           %true time constant, sec
nTicks=round(1000*T/dt);  %ticks ##

% alpha-function/exp filter for past spike counts
range=(10*tauS):-0.1:0;   %past spikes to keep for aplha-filter
% alpha=(range/tauS).*exp(-range/tauS);
alpha=exp(-range/tauS);

%equilibrium potential defined by mean inputs should be > Vthr
vr=(VR+sum(alpha)*dttrue*(scaleE*poisE*VE+scaleI*poisI*VI))/...
  (1+sum(alpha)*dttrue*(scaleE*poisE+scaleI*poisI));
fprintf('Equilibrium state V=%g\n',vr);

%% calc
tic
V=repmat(vr,N,1);             %current membrane voltages
bufE=zeros(N,length(alpha));  %buffer of exc spikes for alpha-filter
bufI=zeros(N,length(alpha));  %buffer of inh spikes for alpha-filter
n=false(N,nTicks);            %raster of spikes
itrack=1;
track=zeros(1,nTicks);        %track voltage of one neuron, for ref
nspk=0;                       %spike counter
fprintf('computing... %8g rate %8g',0,0);
for t=1:nTicks
  %integrated eqn is dV/dt = -A V + B,
  %where A=(1+gE+gI)/tauR, B=(V+gE*VE+gI*VI)/tauR
  %and conductances are alpha-function transformed spike trains
  
  %determine input # of exc spikes for each neuron
  wE=scaleE*(sum(W(:,V>Vthr & flags(:)),2)+poissrnd(poisE*dttrue,N,1));
  
  %determine input # of inh spikes for each neuron
  wI=scaleI*(sum(W(:,V>Vthr & ~flags(:)),2)+poissrnd(poisI*dttrue,N,1));
  
  %update buffers of spike counts
  bufE=[bufE(:,2:end),wE];
  bufI=[bufI(:,2:end),wI];
  
  %determine gE, gI as buf*alpha
  gE=bufE*alpha';
  gI=bufI*alpha';
  
  %constant A & B
  A=(1+gE+gI)/tauR;
  B=(VR+gE*VE+gI*VI)/tauR;
    
  %reset spiked neurons
  nspk=nspk+sum(V>Vthr);
  n(V>Vthr,t)=1;
  V(V>Vthr)=Vres;
  
  %update all neurons
  V=V+dt*(-A.*V+B);
  
  %collect V-track for reference
  track(t)=V(itrack);
  
  if(mod(t,1000)==0) 
    s=repmat('\b',1,22); 
    fprintf([s,'%8g rate %8g'],t*dttrue,nspk/N/(t*dttrue));
  end
end
fprintf('\n');
toc


%% COURTESY PART
if(aim)
  xtime=(1:nTicks)*dttrue;
  figure,plot(xtime,track)
  fprintf('voltage mean/var/rate %g/%g/%g\n',mean(track),std(track),sum(n(itrack,:))/T);
else
  figure,hold on
  for k=1:N 
    idx=find(n(k,:));
    if(flags(k)) s='.'; else s='.r'; end
    plot(idx,repmat(k,size(idx)),s); 
  end
  hold off
end
fprintf('mean rate over pop %g\n',sum(n(:))/N/T);

save result n W flags dt
