function PF_Pop_Inf_v2_5b(T,pl)
clear; %close all; 
T = 10000; pl = 1;
% this script simulates V.Nc cells, and then infers spikes for each,
% assuming they are independent.
%
% 1) set simulation metadata
% 2) initialize parameters
% 3) simulate data
% 4) plot simulation (ie, truth)
% 5) estimate connection matrix from spikes
% 6) prep for spike inference and matrix learning
% 7) omega inference using fast filter
% 8) omega inference using particle filter
% 9) make pretty fig
%
% Remarks:
% a) some of the code is general for V.M spike history terms per neurons, but not all (eg, time constants)
% b) inference assumes all correct parameters (except those governing GLM)
% c) # of external stimulus dimensions is currently restricted to 1
%
% Version Updates:
% v2_1: put EM in loop, assume that each neuron has 1 spike history term in
% E-step, don't plot inference between each EM iteration
%
% v2_2: code is general enough for arbitrary # neurons.  also started
% writing code to estimate matrix using fast method
%
% v2_3: is now a function that takes input T for time. also only saves
% stuff needed for plotting, not other crap.
%
% v2_4: now takes another input, 'pl', which, if 1, plots
% 
% v2_5: renamed simulation structure 'D' from 'S', changed some stuff to
% reduce memory requirements

%% 1) set simulation metadata

pop_sim
V.Nc=V.Ncells;



%% 4) plot simulation results

% if pl==1
%     figure(3), clf, imagesc(P.omega), colormap(gray)
%     tt  = 500;
%     col = [1 0 0; 0 0 1; 0 0.5 0; 1 0.5 0; 1 0 1];          % define colors for mean
%     figure(1), clf
%     for i=1%:V.Nc
%         h1=subplot(311); hold on, plot(((Cell{i}.F(tt:2*tt))./max(Cell{i}.F(tt:2*tt)))+1,'Color',col(i,:)); stem(Cell{i}.n(tt:2*tt),'Color',col(i,:)), axis('tight'), ylabel('F')
%         h2=subplot(312); hold on, plot(y(i,(tt:2*tt)),'Color',col(i,:)), axis('tight'), ylabel('y') %stem(Cell{i}.n,'Color',col(i,:)),
%         h3=subplot(313); hold on, plot(p(i,(tt:2*tt)),'Color',col(i,:)), axis('tight'), ylabel('p')
%     end
%     linkaxes([h1 h2 h3],'x')
%     legend('1','2')
%     for i=1:V.Nc, disp(sum(Cell{i}.n)); end
% end

%% 5) estimate connection matrix directly from spikes

% V.M       = 1;                                    % # spike history terms per neuron (fixed at one for this version of code)
% V.n_params= 1;                                    % if 1, estimate k
% V.h_params= 1;                                    % if 1, estimate omega (self-coupling)
% V.F_params= 0;                                    % if 1, estimate observation parameters
% V.C_params= 0;                                    % whether to compute
% V.StimDim = V.Nc;                               % set external stim dimesions to # cells

V.est_n=1;
V.est_h=1;
V.est_F=0;
V.est_c=0;
V.Nspikehist=1;
V.Nparticles=1;
V.smc_plot=0;
V.smc_iter_max=1;


P.lik=0;
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
    Enew2{i}    = smc_oopsi_m_step(Tim,Cell{i},0,E,Cell{i}.F);
%     Enew2{i}    = GOOPSI_Mstep_v1_0(Tim,Cell{i},0,E,Cell{i}.F);
end

Phat{1}.omega = GetMatrix(V.Nc,Enew2);

if pl==1
    figure(3), clf,
    clims(1)=min(min(P.omega(:)),min(Phat{1}.omega(:)));
    clims(2)=max(max(P.omega(:)),max(Phat{1}.omega(:)));
    subplot(121), imagesc(P.omega,clims), colormap(gray), %colorbar
    subplot(122), imagesc(Phat{1}.omega,clims), %colorbar
    Fs=1024; ts=0:1/Fs:1; sound(sin(2*pi*ts*200)),
    %keyboard
end

%% 6) pop pf preparation

% metaparameters necessary to run smc-em code
% V.N       = 99;                                   % # of particles
V.Mstep   = 1;                                    % whether to estimate parameters
V.pf      = 1;                                    % if 1, then conditional sampler, if 0, then prior sampler
V.freq    = 1;                                    % # time steps per observation (d in BJ08)
V.T_o     = V.T/V.freq;                       % # of observable time steps
V.M       = 1;                                    % # of spike history terms per neuron
V.Mstep   = 0;                                    % whether to estimate parameters
V.ptiles  = 1;                                    % generate percentiles

V.Nparticles=50;
V.fast_do=0;
V.smc_do=1;
% V.StimDim=V.Ncells;

% initialize stuff for each cell
Tim         = V;                                  % copy V structure for input to Mstep function
E           = P;
E.omega     = E.omega(i,i);                         % initialize self-coupling term
E.k         = E.k*ones(V.StimDim,1);              % initialize external stim and cross-coupling terms
for i=1:V.Nc,
    I{i}.M.nbar = zeros(1,V.T);                   % initialize spike trains
    I{i}.P      = E;                                % initialize parameters
    smc{i}.P    = E;
end

%% 7) Pop inference using fast filter

% E.tau       = E.tau_c;                                  % calcium decay time constant (sec)
% E.sig       = E.sigma_c;                                    % noise on
% E.lam       = V.T/(rate*E.A)*V.dt;              % expected jump size ber time bin
% V.Plot    = 1;
% V.MaxIter = 10;                                   % max number of EM iterations
%
% for i=1:V.Nc,                                     % initialize estimate of spike trains
%     [Fn{i} FP{i}]       = FOOPSI_v1_9(Cell{i}.F',E,V);
% end
%
% for i=1:V.Nc
%
%     % append external stimulus for neuron 'i' with spike histories from other cells
%     h = zeros(V.Nc-1,V.T);                  % we append this to x to generate input into neuron from other neurons
%     Pre=1:V.Nc;                               % generate list of presynaptic neurons
%     Pre(Pre==i)=[];                             % remove self
%     k=0;                                        % counter of dimension
%     for j=Pre                                   % loop thru all presynaptic neurons
%         k=k+1;                                  % generate input to neuron based on posterior mean spike train from neuron j
%         h(k,:) = filter(1,[1 -(1-V.dt/P.tau_h)],Fn{j}/max(Fn{j}));
%     end
%     Tim.x   = [V.x; h];                         % append input from other neurons onto external stimulus
%     Tim.N   = 1;
%     T.w_b   = ones(1,V.T);
%     T.n     = Fn{i};
%     T.C     = filter(1,[1 -(1-V.dt/FP{i}.tau)],Fn{j}/max(Fn{j}));
%     T.h     = filter(1,[1 -(1-V.dt/E.tau_h)],Fn{j}/max(Fn{j}));;
%     FP{i}   = GOOPSI_Mstep_v1_0(Tim,T,0,II{i}.P,Cell{i}.F);
% end
% PlotPop_v1_0(V,FE,P.omega,Phat{1}.omega)

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
        smc{i}  = run_oopsi(Cell{i}.F,V,smc{i}.P);
%         [S M II{i}.P]   = GOOPSI_main_v3_0(Cell{i}.F,I{i}.P,Tim);
        II{i}.S.n       = smc{i}.E.n;
        II{i}.S.h       = smc{i}.E.h; 
        II{i}.S.w_b     = smc{i}.E.w; 
        II{i}.M.nbar    = smc{i}.E.nbar;
        II{i}.P         = smc{i}.P;
        smc{i}.E.w_b    = smc{i}.E.w;
%         II{i}.S.n       = S.n; 
%         II{i}.S.h       = S.h; 
%         II{i}.S.w_b     = S.w_b; 
%         II{i}.M.nbar    = M.nbar;
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
        II{i}.P.k=smc{i}.P.k*ones(V.Ncells,1);
        smc{i}.P = smc_oopsi_m_step(Tim,smc{i}.E,0,II{i}.P,Cell{i}.F);
        EE{i}    = smc{i}.P;
    end
    Phat{2}.omega = GetMatrix(V.Nc,EE);
    PlotMatrix(Cell,V,P,Phat)
    save(['VConnector', num2str(V.T), '_', num2str(V.Nc), '_', num2str(tr)],'Phat','D','V','P')
    Fs=tr*1024; ts=0:1/Fs:1; sound(sin(2*pi*ts*200)),
end

end

function omega = GetMatrix(Nc,P)

omega=zeros(Nc);
for i=1:Nc
    omega(i,i)=P{i}.omega;
    Pre=1:Nc;                                       % generate list of presynaptic neurons
    Pre(Pre==i)=[];                                 % remove self
    k=0;                                            % counter of dimension
    for j=Pre
        k=k+1;                                      % generate input to neuron based on posterior mean spike train from neuron j
        omega(i,j)=P{i}.k(k+1);
    end
end

end

function PlotMatrix(Cell,V,P,Phat)

%%
fig=figure(5); clf,
col = [1 0 0; 0 0 1; 0 0.5 0; 1 0.5 0; 1 0 1];          % define colors for mean

% Phat.omega=zeros(V.Nc);
% for i=1:V.Nc
%     Phat.omega(i,i)=I{i}.P.omega;
%     Pre=1:V.Nc;                                   % generate list of presynaptic neurons
%     Pre(Pre==i)=[];                                 % remove self
%     k=0;                                            % counter of dimension
%     for j=Pre
%         k=k+1;                                      % generate input to neuron based on posterior mean spike train from neuron j
%         Phat.omega(i,j)=I{i}.P.k(k+1);
%     end
% end

tstart=2/V.dt;
tdt=1/V.dt;
tend = min(tstart+8*tdt,V.T);
subplot(2,3,[1 2 3]), hold on
for i=1:2
    plot(((Cell{i}.F(tstart:tend))./(Cell{i}.F(tstart:tend)))+1,'Color',col(i,:));
end
axis('tight'),
xticks  = tstart:tdt:tend;               % XTick positions
set(gca,'YTick',[],'YTickLabel',[])
set(gca,'XTick',xticks,'XTickLabel',xticks*V.dt)
ylabel('Fluorescence')
xlabel('Time (sec)')

clims(1)=min([P.omega(:)' Phat{1}.omega(:)' Phat{2}.omega(:)']);
clims(2)=max([P.omega(:)' Phat{1}.omega(:)' Phat{2}.omega(:)']);
subplot(234), imagesc(P.omega,clims), colormap(gray),
set(gca,'XTick',[1:V.Nc],'YTick',[1:V.Nc]) %colorbar
title('True matrix')
ylabel('Presynaptic'), xlabel('Postsynaptic')

subplot(235), imagesc(Phat{1}.omega,clims), %colorbar
set(gca,'XTick',[1:V.Nc],'YTick',[1:V.Nc]) %colorbar
title('Matrix from spikes')
% ylabel('Presynaptic'), xlabel('Postsynaptic')

subplot(236), imagesc(Phat{2}.omega,clims), %colorbar
set(gca,'XTick',[1:V.Nc],'YTick',[1:V.Nc]) %colorbar
title('Matrix from fluorescence')
% ylabel('Presynaptic'), xlabel('Postsynaptic')

% print fig
wh=[7 5];   %width and height
set(fig,'PaperPosition',[0 11-wh(2) wh]);
print('-depsc','VConnector')
print('-dpdf','VConnector')

end