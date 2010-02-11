
%% simulat or load data

pl          = 1;
pop_sim

%% 5) estimate connection matrix directly from spikes

V.M       = 1;                                    % # spike history terms per neuron (fixed at one for this version of code)
V.n_params= 1;                                    % if 1, estimate k
V.h_params= 1;                                    % if 1, estimate omega (self-coupling)
V.F_params= 0;                                    % if 1, estimate observation parameters
V.C_params= 0;                                    % whether to compute
V.StimDim = V.Nc;                               % set external stim dimesions to # cells
Tim         = V;                                  % Tim is V for this estimation
Tim.N       = 1;                                    % # of particles
for i=1:V.Nc
    h=zeros(V.Nc-1,V.T);
    Pre=1:V.Nc;                                   % generate list of presynaptic neurons
    Pre(Pre==i)=[];                                 % remove self
    k=0;                                            % counter of dimension
    for j=Pre                                       % loop thru all presynaptic neurons
        k=k+1;                                      % generate input to neuron based on posterior mean spike train from neuron j
        h(k,:) = Cell{j}.h;
    end
    Tim.x       = [V.x; h];                       % append input from other neurons onto external stimulus
    E           = P;
    E.omega     = E.omega(i,i);                     % initialize self-coupling term
    E.k         = E.k*ones(V.Nc,1);               % initialize external stim and cross-coupling terms
    Cell{i}.w_b    = ones(1,V.T);
    Enew2{i} = smc_oopsi_m_step(Tim,Cell{i},0,E,Cell{i}.F);
%     Enew2{i}    = GOOPSI_Mstep_v1_0(Tim,Cell{i},0,E,Cell{i}.F);
end

Phat{1}.omega = GetMatrix2(V.Nc,Enew2);

% if pl==1
%     figure(3), clf,
%     clims(1)=min(min(P.omega(:)),min(Phat{1}.omega(:)));
%     clims(2)=max(max(P.omega(:)),max(Phat{1}.omega(:)));
%     subplot(121), imagesc(P.omega,clims), colormap(gray), %colorbar
%     subplot(122), imagesc(Phat{1}.omega,clims), %colorbar
%     Fs=1024; ts=0:1/Fs:1; sound(sin(2*pi*ts*200)),
% end

%% 6) pop pf preparation

% metaparameters necessary to run smc-em code
V.N       = 99;                                   % # of particles
V.Mstep   = 1;                                    % whether to estimate parameters
V.pf      = 1;                                    % if 1, then conditional sampler, if 0, then prior sampler
V.freq    = 1;                                    % # time steps per observation (d in BJ08)
V.T_o     = V.T/V.freq;                       % # of observable time steps
V.M       = 1;                                    % # of spike history terms per neuron
V.Mstep   = 0;                                    % whether to estimate parameters
V.ptiles  = 1;                                    % generate percentiles

% initialize stuff for each cell
Tim         = V;                                  % copy V structure for input to Mstep function
E           = P;
E.omega     = E.omega(i,i);                         % initialize self-coupling term
E.k         = E.k*ones(V.StimDim,1);              % initialize external stim and cross-coupling terms
for i=1:V.Nc,
    I{i}.M.nbar = zeros(1,V.T);                   % initialize spike trains
    I{i}.P      = E;                                % initialize parameters
end

%% 8)  Pop EM iteration

for tr=1:1
    fprintf('\nTrial # %g\n',tr)

    % for each neuron, infer spike train conditioned on previous EM
    % iterations spike history terms
    for i=1:V.Nc,                                 % infer spikes for each neuron

        % append external stimulus for neuron 'i' with spike histories from other cells
        h = zeros(V.Nc-1,V.T);                  % we append this to x to generate input into neuron from other neurons
        Pre=1:V.Nc;                               % generate list of presynaptic neurons
        Pre(Pre==i)=[];                             % remove self
        k=0;                                        % counter of dimension
        for j=Pre                                   % loop thru all presynaptic neurons
            k=k+1;                                  % generate input to neuron based on posterior mean spike train from neuron j
            h(k,:) = filter(1,[1 -(1-V.dt/P.tau_h)],I{j}.M.nbar);
        end
        Tim.x = [V.x; h];                         % append input from other neurons onto external stimulus

        % infer spike train for neuron 'i'
        fprintf('\nNeuron # %g\n',i)
        if 0
            [S M II{i}.P]   = fast_oopsi(Cell{i}.F,Tim,I{i}.P);
            II{i}.S.n       = S ./ max(S)';
            II{i}.S.h       = zeros(1,V.T);
            for t=2:Tim.T
                II{i}.S.h(t)=(1-V.dt/P.tau_h)*II{i}.S.h(t-1)+II{i}.S.n(t-1);
            end
            II{i}.S.w_b     = ones(1,V.T);
            II{i}.M.nbar    = S ./ max(S)';
            II{i}.P.k       = P.k*ones(V.Nc,1);
            II{i}.P.omega   = 0;
            Tim.N = 1;
        else
            [S M II{i}.P]   = GOOPSI_main_v3_0(Cell{i}.F,I{i}.P,Tim);
            II{i}.S.n       = S.n;
            II{i}.S.h       = S.h;
            II{i}.S.w_b     = S.w_b;
            II{i}.M.nbar    = M.nbar;
        end
    end

    % set inference for each neuron to the newly updated inference
    for i=1:V.Nc, I{i}.S = II{i}.S; I{i}.M = II{i}.M; end

    % given new inference for each neuron, update parameters
    for i=1:V.Nc

        % append external stimulus for neuron 'i' with spike histories from other cells
        h = zeros(V.Nc-1,V.T);                  % we append this to x to generate input into neuron from other neurons
        Pre=1:V.Nc;                               % generate list of presynaptic neurons
        Pre(Pre==i)=[];                             % remove self
        k=0;                                        % counter of dimension
        for j=Pre                                   % loop thru all presynaptic neurons
            k=k+1;                                  % generate input to neuron based on posterior mean spike train from neuron j
            h(k,:) = filter(1,[1 -(1-V.dt/P.tau_h)],I{j}.M.nbar);
        end
        Tim.x = [V.x; h];                         % append input from other neurons onto external stimulus

        fprintf('\nNeuron # %g\n',i)
        I{i}.P = smc_oopsi_m_step(Tim,I{i}.S,0,II{i}.P,Cell{i}.F);
%         I{i}.P  = GOOPSI_Mstep_v1_0(Tim,I{i}.S,0,II{i}.P,Cell{i}.F);
        EE{i}    = I{i}.P;
    end
    Phat{2}.omega = GetMatrix2(V.Nc,EE);
    PlotMatrix2(Cell,V,P,Phat)
    save(['SimConnector', num2str(V.T), '_', num2str(V.Nc), '_', num2str(tr)],'Phat','Cell','V','P')
    Fs=tr*1024; ts=0:1/Fs:1; sound(sin(2*pi*ts*200)),
end


