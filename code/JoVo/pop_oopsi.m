clear, clc
%% simulat or load data
% (see pop_sim to know which fields of V must be set here)

V.name='sim1';

pop_sim

%% set up variables

V.est_n=1;          % estimate spiking parameters
V.est_h=1;          % estimate spike history parameters
V.est_F=0;          % do NOT estimate fluorescence parameters
V.est_c=0;          % do NOT estimate calcium parameters
V.Nspikehist=1;     % assume 1 spike history terms
V.Nparticles=50;    % # of monte carlo samples
V.smc_plot=0;       % do not plot smc
V.smc_iter_max=1;   % number of iterations
V.fast_do=0;        % do not do fast-oopsi
V.smc_do=1;         % do not do smc-oopsi

%% estimate connection matrix directly from spikes

use_spikes=1;
if use_spikes
    Phat{1}.omega = GetWSpikes(Cell,V,P);
end

%% pop pf preparation

% initialize stuff for each cell
V.StimDim   = V.Ncells;                 % set external stim dimesions to # cells
Tim         = V;                        % copy V structure for input to Mstep function
E           = P;
E.omega     = E.omega(i,i);             % initialize self-coupling term
E.k         = E.k*ones(V.StimDim,1);    % initialize external stim and cross-coupling terms
for i=1:V.Ncells,
    I{i}.M.nbar = zeros(1,V.T);         % initialize spike trains
    I{i}.P      = E;                    % initialize parameters
    smc{i}.P    = E;
end

for tr=1:1                                  % iterate trials

    % infer spikes for each neuron
    for i=1:V.Ncells,                       
        
        % append external stimulus for neuron 'i' with spike histories from other cells
        h = zeros(V.Ncells-1,V.T);          % we append this to x to generate input into neuron from other neurons
        Pre=1:V.Ncells;                     % generate list of presynaptic neurons
        Pre(Pre==i)=[];                     % remove self
        k=0;                                % counter of dimension
        for j=Pre                           % loop thru all presynaptic neurons
            k=k+1;                          % generate input to neuron based on posterior mean spike train from neuron j
            h(k,:) = filter(1,[1 -(1-V.dt/P.tau_h)],I{j}.M.nbar);
        end
        Tim.x = [V.x; h];                   % append input from other neurons onto external stimulus

        % infer spike train for neuron 'i'
        fprintf('\nNeuron # %g\n',i)
        smc{i}          = run_oopsi(Cell{i}.F,Tim,smc{i}.P);
        II{i}.P         = smc{i}.P;
        II{i}.S.n       = smc{i}.E.n;
        II{i}.S.h       = smc{i}.E.h;
        II{i}.S.w_b     = smc{i}.E.w;
        II{i}.M.nbar    = smc{i}.E.nbar;
    end

    % set inference for each neuron to the newly updated inference
    for i=1:V.Ncells,
        I{i}.S = II{i}.S;
        I{i}.M = II{i}.M;
    end

    % estimate connectivity
    for i=1:V.Ncells

        % append external stimulus for neuron 'i' with spike histories from other cells
        h = zeros(V.Ncells-1,V.T);          % we append this to x to generate input into neuron from other neurons
        Pre=1:V.Ncells;                     % generate list of presynaptic neurons
        Pre(Pre==i)=[];                     % remove self
        k=0;                                % counter of dimension
        for j=Pre                           % loop thru all presynaptic neurons
            k=k+1;                          % generate input to neuron based on posterior mean spike train from neuron j
            h(k,:) = filter(1,[1 -(1-V.dt/P.tau_h)],I{j}.M.nbar);
        end
        Tim.x = [V.x; h];                   % append input from other neurons onto external stimulus

        fprintf('\nNeuron # %g\n',i)
        I{i}.P = smc_oopsi_m_step(Tim,I{i}.S,0,II{i}.P,Cell{i}.F);
        EE{i}    = I{i}.P;
    end
    Phat{tr+1}.omega = GetMatrix(V.Ncells,EE);
end
PlotMatrix(Cell,V,P,Phat)
save(['../../data/' V.name] ,'Phat','Cell','V','P')
Fs=1024; ts=0:1/Fs:1; sound(sin(2*pi*ts*200)),