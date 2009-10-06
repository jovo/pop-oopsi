function S = gIAF_fast2
% this function takes no input, but calls gIAF_params
% (which is appended below) to determine the cell parameters
% the nomenclature is used to be in accordance with JHDC, such that
% this function, or these parameters can be simply copied and pasted
% directly into a JHDC worker function.  the neuron model is a
% conductance-based integrate-and-fire.  we use a propogator matrix to
% calculate the conductance at every step, and use a numerical exponential
% euler to solve for voltage

clear, clc
starttime=cputime;  %get start time to keep track of how long simulation took
P = gIAF_params;    %set parameters
MaxSpikeNum = P.nSec/P.dt;


S       = [];       %Extize the data I want to save in a structure
S.t     = 1;        %keep track of time
S.dt    = P.dt;     %step size

S.E.V                   = zeros(5,P.E_Ncells);  %voltage Extization
S.E.V(5,:)              = -rand(P.E_Ncells,1)*abs(mean(P.E_Vth')-P.E_Vreset)+mean(P.E_Vth);               %set last component to Vm
S.E.last                = P.E_last;             %list of the most recent spike times
S.E.SpikeTimes          = zeros(1,MaxSpikeNum);         %list of spike times
S.E.SpikeCellNums       = zeros(1,MaxSpikeNum);         %list of which cell spiked at each time
S.E.nSpikes             = 0;                    %cumulative # of excitatory spikes
S.E.spiked{P.E_Delay}   = [];                   %cells that spiked in this time step

S.I.V                   = zeros(5,P.I_Ncells);
S.I.V(5,:)              = P.I_Vm;
S.I.last                = P.I_last;
S.I.SpikeTimes          = zeros(1,1e3);
S.I.SpikeCellNums       = zeros(1,1e3);
S.I.nSpikes             = 0;
S.I.spiked{P.I_Delay}   = [];

S.E.Varch               = zeros(P.E_Ncells, P.nSec/S.dt);   %archive voltage for each cell
S.I.Varch               = zeros(P.I_Ncells, P.nSec/S.dt);   %archive voltage for each cell
% S.Vrand                 = .2*randn(P.Ncells, P.nSec/S.dt);   %archive voltage for each cell

S.E.Ge_arch               = zeros(P.E_Ncells, P.nSec/S.dt);   %archive voltage for each cell
S.I.Ge_arch               = zeros(P.I_Ncells, P.nSec/S.dt);   %archive voltage for each cell
S.E.Gi_arch               = zeros(P.E_Ncells, P.nSec/S.dt);   %archive voltage for each cell
S.I.Gi_arch               = zeros(P.I_Ncells, P.nSec/S.dt);   %archive voltage for each cell

% determine weights for all the cells
% load PreNestWeights
P.EPSP.scale=0.05;
P.IPSP.scale=15;
[W,s]=make_weights(0.1,P);
W=W*0.01;
N=P.Ncells;
Nexc=P.E_Ncells;
Ninh=P.I_Ncells;

for n=1:Nexc                              %this loop determines all the neurons that are postsynaptic to the excitatory population
    EE.PrePost{n} = find(W(1:Nexc,n));
    EE.Weights{n} = 0*W(EE.PrePost{n},n);

    IE.PrePost{n} = find(W(Nexc+1:N,n));  %similar, but stores all the post-synaptic neurons to the inhibitory population
    IE.Weights{n} = 0*W(Nexc+IE.PrePost{n},n);
end

for n=1:Ninh                              %loop through all the neurons that are presynaptic to the inhibitory population
    EI.PrePost{n} = find(W(1:Nexc,Nexc+n));        %E's that are post to the I's
    EI.Weights{n} = 0*5*W(EI.PrePost{n},Nexc+n);

    II.PrePost{n} = find(W(Nexc+1:N,Nexc+n));      %I's that are post to the I's
    II.Weights{n} = 0*0*W(Nexc+II.PrePost{n},Nexc+n);
end

while S.t<(P.nSec/S.dt)

    %Force some neurons to spike
%     if any(S.t==P.Ext.ForcedSpikeTimes)
%         ForcedSpikes = P.Ext.ForcedCellNums(find(S.t==P.Ext.ForcedSpikeTimes));
%     else
        ForcedSpikes = [];
%     end
    
%     if S.t>P.Ext.Init.Start & S.t<=(P.Ext.Init.Start + P.Ext.Init.Duration)...
%             & mod(S.t,round(1/(P.Ext.Init.Freq*S.dt)))==0
%         ForcedSpikes = P.Ext.Init.CellNums; %store a vector which corresponds to the neurons that were forced to spike in this time step
%     else
%         ForcedSpikes = [];
%     end

    E.Ge = P.E_Pe * S.E.V(1:2,:);
    E.Gi = P.E_Pi * S.E.V(3:4,:);

    I.Ge = P.I_Pe * S.I.V(1:2,:);
    I.Gi = P.I_Pi * S.I.V(3:4,:);
   
    S.E.Ge_arch(:,S.t) = E.Ge(2,:);
    S.E.Gi_arch(:,S.t) = E.Gi(2,:);
    
    S.I.Ge_arch(:,S.t) = I.Ge(2,:);
    S.I.Gi_arch(:,S.t) = I.Gi(2,:);

    SEspiked=S.E.spiked{1};
    for n=1:length(SEspiked) %add alpha function of size determined by weights to each of the E voltages
        S.E.V(1,EE.PrePost{SEspiked(n)}) = S.E.V(1,EE.PrePost{SEspiked(n)}) + EE.Weights{SEspiked(n)}'; %add input to 1st component of vector
        S.I.V(1,IE.PrePost{SEspiked(n)}) = S.I.V(1,IE.PrePost{SEspiked(n)}) + IE.Weights{SEspiked(n)}';
    end

    SIspiked=S.I.spiked{1};
    for m=1:length(SIspiked)
        S.E.V(3,EI.PrePost{SIspiked(m)}) = S.E.V(3,EI.PrePost{SIspiked(m)}) + EI.Weights{SIspiked(m)}';
        S.I.V(3,II.PrePost{SIspiked(m)}) = S.I.V(3,II.PrePost{SIspiked(m)}) + II.Weights{SIspiked(m)}';
    end

    S.E.V(5,:) = Integrate(S.E.V(5,:), E, P, 'E');
    S.I.V(5,:) = Integrate(S.I.V(5,:), I, P, 'I');

    S.E.V(5,:) = S.E.V(5,:);% + S.Vrand(1:P.E_Ncells,S.t)';
    S.I.V(5,:) = S.I.V(5,:);% + S.Vrand(P.E_Ncells+1:P.Ncells,S.t)';
    
    S.E.Varch(:,S.t)=S.E.V(5,:);
    S.I.Varch(:,S.t)=S.I.V(5,:);

    S.E = Fire(S, P, 'E', ForcedSpikes); %determine which cells are refractory, which have spiked, and adjust voltages accordingly
    S.I = Fire(S, P, 'I');

    S.t=S.t+1;  %step

    for n=1:P.E_Delay-1
        S.E.spiked{n}=S.E.spiked{n+1};
    end

    for n=1:P.I_Delay-1
        S.I.spiked{n}=S.I.spiked{n+1};
    end
end

plot_hans(S,P,2)    %plot data
S.time=cputime-starttime;
S.W.EE=EE;
S.w=W;

end

%% this function generates all the parameters
function P = gIAF_params

% general parametes
P.dt        = 1;  %time step size
P.nSec      = 5000; %duration of simulation in ms
P.Ncells    = 50;  %total number of cells
P.Nrepeats  = 2;

% PSP properties
P.EPSP.tauSyn  = 1;  %synaptic time constant
P.IPSP.tauSyn  = 2;

% Excitatory cell parameters
P.E_Ncells  = 0.8 * P.Ncells; %# of excitatory cells
P.E_Vm      = -35;            %resting membrane potential must be zero for the math to work
P.E_Vth     = rand(1,P.E_Ncells)*(-10)-35;          %threshold
P.E_Vreset  = -65;       %reset potential
P.E_Ref     = 5;            %absolute refractory period
P.E_last    = zeros(P.E_Ncells,1)-1000; %P.E_Ref/P.dt-1;  %list of the most recent spike times
P.E_tauMem  = 100;            %membrane time constant
P.E_gm      = 1/100;         %conductance
P.E_Rev     = 0;
P.E_Delay   = 50;
[P.E_Pe P.E_Pi] = AddPropagator(P.EPSP.tauSyn, P.IPSP.tauSyn, P.E_tauMem, P.E_gm, P.dt);

% Inhibitory cell parameters
P.I_Ncells  = P.Ncells-P.E_Ncells;
P.I_Vm     = -50;
P.I_Vth    = -45;
P.I_Vreset = -65;
P.I_Ref     = 1;
P.I_last    = zeros(P.I_Ncells,1)-1000;%P.I_Ref/P.dt-1; %10*[1:P.I_Ncells];%
P.I_tauMem  = 100;
P.I_gm      = 1/100;
P.I_Rev     = -90;
P.I_Delay   = 5;
[P.I_Pe P.I_Pi] = AddPropagator(P.EPSP.tauSyn, P.IPSP.tauSyn, P.I_tauMem, P.I_gm, P.dt);

%weight matrix parameters
P.ScaleIE        = 6e5 / P.E_Ncells / P.Nrepeats;  %scale weights from X to Y 
P.ScaleII        = 0   / P.I_Ncells / P.Nrepeats;  %scale weights from Y to Y
P.ScaleEI        = 1e3 / P.I_Ncells / P.Nrepeats;  %scale weights from Y to X

%initial conditions
P.Ext.Init.Start    = 0;        %time to start spike train
P.Ext.Init.Duration = 5/P.dt;   %time to stop spike train
P.Ext.Init.Freq     = 10;       %freq of spike train (in ms^-1)
P.Ext.Init.CellNums = 1;        %cell #'s to spike train

% P.Ext.ForcedSpikeTimes  = [0.1 10 20 30 40] / P.dt;
% P.Ext.ForcedCellNums    = 1*ones(size(P.Ext.ForcedSpikeTimes));

P.Ext.ForcedSpikeTimes  = [];
P.Ext.ForcedCellNums    = [];
for k=1:5
  ind=find(rand(1,P.nSec/P.dt)<0.005*P.dt);
  P.Ext.ForcedSpikeTimes  = [P.Ext.ForcedSpikeTimes ind];
  P.Ext.ForcedCellNums    = [P.Ext.ForcedCellNums repmat(k,1,length(ind))];
end

end

%% this function generates the propogator matrix
function [Pe Pi] = AddPropagator(tauA, tauB, tauM, gm, dt)

PA = exp(-dt/tauA);
PB = exp(-dt/tauB);
PM  = exp(-dt/tauM);

c=1/(tauM*gm);

Pe = [PA, 0; dt*PA, PA];
Pi = [PB, 0; dt*PB, PB];

end

%% this function integrates for a conductance based IAF using exponential euler's method
function V = Integrate(V, G, P, EI)

if EI == 'E'
    Vm = P.E_Vm;
    Gm = P.E_gm;
    Tau = P.E_tauMem;
    Ncells = P.E_Ncells;
elseif EI == 'I'
    Vm = P.I_Vm;
    Gm = P.I_gm;
    Tau = P.I_tauMem;
    Ncells = P.I_Ncells;
else
    error('cell must be E or I, schmuck')
end

Ge=G.Ge(2,:);   %calc excitatory conductance for cells
Gi=G.Gi(2,:);   %calc inhibitory conductance for cells

A=(Vm+(Ge*P.E_Rev+Gi*P.I_Rev)/Gm)/Tau;  %exponential euler's method for integration
B=(1+(Ge+Gi)/Gm)/Tau;
Bexp=exp(-B*P.dt);                      %makes computation slightly faster

V=V.*Bexp+(A./B).*(1-Bexp);             %calc voltage for x cells at next time step

end

%% this function tests for refractoriness and spikes
function SO = Fire(S, P, EI, ForcedSpikes)

if EI == 'E'
    SO      = S.E;
    Ref     = P.E_Ref;
    Vreset  = P.E_Vreset;
    Vth     = P.E_Vth;
    Delay   = P.E_Delay;
elseif EI == 'I'
    SO      = S.I;
    Ref     = P.I_Ref;
    Vreset  = P.I_Vreset;
    Vth     = P.I_Vth;
    Delay   = P.I_Delay;
else
    error('Cell must be E or I, schmuck')
end

%check for refractoriness
SO.ref          = S.t-SO.last<(Ref/S.dt);
SO.V(5,SO.ref)  = Vreset; %on the last component of vector for the V

%spike
if nargin == 4  %since i only force excitatory spikes, i only input ForceSpikes for excitatory cells
    SO.spiked{Delay} = [find(SO.V(5,:)>Vth), ForcedSpikes];
else
    SO.spiked{Delay} = find(SO.V(5,:)>Vth);
end
SO.V(5,SO.spiked{Delay}) = Vreset;                                     %set voltage to reset value
len = length(SO.spiked{Delay});                                        %compute length only once
SO.SpikeCellNums(SO.nSpikes+1 : SO.nSpikes+len) = SO.spiked{Delay};    %store the cell's that spiked
SO.SpikeTimes(SO.nSpikes+1 : SO.nSpikes+len) = S.t;             %store the time for the spikes
SO.nSpikes = SO.nSpikes + len;                                  %increment the scalar cumulative spikes;
SO.last(SO.spiked{Delay}) = S.t;                                       %update the last spike for each neuron

end

%% this function plots the data
function X = plot_hans(S,P,a)
%% the input all the parameters (S and P) and the figure # (a)

figure(a), clf, hold on
plot(S.E.SpikeTimes*S.dt, S.E.SpikeCellNums,'.b')
plot(S.I.SpikeTimes*S.dt, S.I.SpikeCellNums + P.E_Ncells,'.r')
hold off

end

%% this function makes the weights
function [w S] = make_weights(SP,P)

N=P.Ncells;
Ninh=P.I_Ncells;

S=zeros(1,N);
w=sprand(N,N,SP);

%excitatory neurons -- log-rates
idx=1:N-Ninh; S(idx)=1;
W=w(:,idx); ind=find(W>0);
z=P.EPSP.scale*exprnd(1,length(ind),1);
W(ind)=z; w(:,idx)=W;

%inhibitory neurons -- log-rates
idx=N-Ninh+1:N; S(idx)=-1;
W=w(:,idx); ind=find(W>0);
z=P.IPSP.scale*exprnd(1,length(ind),1);
W(ind)=z; w(:,idx)=W;

%remove self-coupling, refractory is Omega, not W
for k=1:N w(k,k)=0; end

end