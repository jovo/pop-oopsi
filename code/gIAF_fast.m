function S = gIAF_fast
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
S.E.V(5,:)              = P.E_Vm;               %set last component to Vm
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

% determine weights for all the cells
load PreNestWeights

for n=1:P.E_Ncells                          %this loop determines all the neurons that are presynaptic to the excitatory population
    EE.PrePost{n} = find(W.EE(:,n));
    EE.Weights{n} = W.EE(EE.PrePost{n},1);

    IE.PrePost{n} = 1:P.I_Ncells;           %similar, but stores all the post-synaptic neurons to the inhibitory population
    IE.Weights{n} = P.ScaleIE;
end

for n=1:P.I_Ncells                          %loop through all the neurons that are presynaptic to the inhibitory population
    EI.PrePost{n} = 1:P.E_Ncells;           %E's that are pre to the I's
    EI.Weights{n} = P.ScaleEI;

    II.PrePost{n} = 1:P.I_Ncells;           %I's that are pre to the E's
    II.Weights{n} = P.ScaleII;
end

while S.t<(P.nSec/S.dt)

    %Force some neurons to spike
    if any(S.t==P.Ext.ForcedSpikeTimes)
        ForcedSpikes = P.Ext.ForcedCellNums(find(S.t==P.Ext.ForcedSpikeTimes));
    else
        ForcedSpikes = [];
    end
    
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

    S.E.Varch(:,S.t)=S.E.V(5,:);

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

end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% this function generates all the parameters
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function P = gIAF_params

% general parametes
P.dt        = 0.1;  %time step size
P.nSec      = 500; %duration of simulation in ms
P.Ncells    = 200;  %total number of cells
P.Nrepeats  = 2;

% PSP properties
P.EPSP.tauSyn  = 2;  %synaptic time constant
P.IPSP.tauSyn  = 2;

% Excitatory cell parameters
P.E_Ncells  = P.Ncells;     %# of excitatory cells
P.E_Vm      = -50;            %resting membrane potential must be zero for the math to work
P.E_Vth     = -45;          %threshold
P.E_Vreset  = -55;       %reset potential
P.E_Ref     = 5;            %absolute refractory period
P.E_last    = zeros(P.E_Ncells,1)-1000; %P.E_Ref/P.dt-1;  %list of the most recent spike times
P.E_tauMem  = 65;            %membrane time constant
P.E_gm      = 1/50;         %conductance
P.E_Rev     = 52.5;
P.E_Delay   = 50;
[P.E_Pe P.E_Pi] = AddPropagator(P.EPSP.tauSyn, P.IPSP.tauSyn, P.E_tauMem, P.E_gm, P.dt);

% Inhibitory cell parameters
P.I_Ncells  = 0.1 * P.Ncells;
P.I_Vm     = -59;
P.I_Vth    = -55;
P.I_Vreset = -65;
P.I_Ref     = 1;
P.I_last    = zeros(P.I_Ncells,1)-1000;%P.I_Ref/P.dt-1; %10*[1:P.I_Ncells];%
P.I_tauMem  = 100;
P.I_gm      = 1/100;
P.I_Rev     = -80;
P.I_Delay   = 5;
[P.I_Pe P.I_Pi] = AddPropagator(P.EPSP.tauSyn, P.IPSP.tauSyn, P.I_tauMem, P.I_gm, P.dt);

%weight matrix parameters
P.ScaleIE        = 0*6e2 / P.E_Ncells / P.Nrepeats;  %scale weights from X to Y 
P.ScaleII        = 0   / P.I_Ncells / P.Nrepeats;  %scale weights from Y to Y
P.ScaleEI        = 0*1e1 / P.I_Ncells / P.Nrepeats;  %scale weights from Y to X

%initial conditions
P.Ext.Init.Start    = 0;        %time to start spike train
P.Ext.Init.Duration = 5/P.dt;   %time to stop spike train
P.Ext.Init.Freq     = 10;       %freq of spike train (in ms^-1)
P.Ext.Init.CellNums = 1;        %cell #'s to spike train

P.Ext.ForcedSpikeTimes  = [0.1 10 20 30 40] / P.dt;
P.Ext.ForcedCellNums    = 1*ones(size(P.Ext.ForcedSpikeTimes));

end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% this function generates the propogator matrix
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Pe Pi] = AddPropagator(tauA, tauB, tauM, gm, dt)

PA = exp(-dt/tauA);
PB = exp(-dt/tauB);
PM  = exp(-dt/tauM);

c=1/(tauM*gm);

Pe = [PA, 0; dt*PA, PA];
Pi = [PB, 0; dt*PB, PB];

end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% this function integrates for a conductance based IAF
% %% using exponential euler's method
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% this function tests for refractoriness and spikes
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% this function plots the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function X = plot_hans(S,P,a)
%% the input all the parameters (S and P) and the figure # (a)

figure(a), clf, hold on
plot(S.E.SpikeTimes*S.dt, S.E.SpikeCellNums,'.b')
plot(S.I.SpikeTimes*S.dt, S.I.SpikeCellNums + P.E_Ncells,'.r')
hold off

end