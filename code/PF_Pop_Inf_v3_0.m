function [S M E I II] = PF_Pop_Inf_v3_0(Sim,P,F)
% this script simulates Sim.Nc cells, and then infers spikes for each,
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
% Input:
%   T:  # of time steps (min of 1000, when plotting)
%   pl: whether to plot stuff between iterations
% 
% Remarks:
% a) some of the code is general for Sim.M spike history terms per neurons, but not all (eg, time constants)
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
% 
% v3_0: takes input Sim,P,F, and outputs 
%%

% initialize stuff for each cell
Tim         = Sim;                                  % copy Sim structure for input to Mstep function
E           = P;
E.omega     = E.omega(i,i);                         % initialize self-coupling term
E.k         = E.k*ones(Sim.StimDim,1);              % initialize external stim and cross-coupling terms
for i=1:Sim.Nc,
    I{i}.M.nbar = zeros(1,Sim.T);                   % initialize spike trains
    I{i}.P      = E;                                % initialize parameters
end

%% loop over inference
for tr=1:5
    fprintf('\nTrial # %g\n',tr)

    % for each neuron, infer spike train conditioned on previous EM
    % iterations spike history terms
    for i=1:Sim.Nc,                                 % infer spikes for each neuron

        % append external stimulus for neuron 'i' with spike histories from other cells
        h = zeros(Sim.Nc-1,Sim.T);                  % we append this to x to generate input into neuron from other neurons
        Pre=1:Sim.Nc;                               % generate list of presynaptic neurons
        Pre(Pre==i)=[];                             % remove self
        k=0;                                        % counter of dimension
        for j=Pre                                   % loop thru all presynaptic neurons
            k=k+1;                                  % generate input to neuron based on posterior mean spike train from neuron j
            h(k,:) = filter(1,[1 -(1-Sim.dt/P.tau_h)],I{j}.M.nbar);
        end
        Tim.x = [Sim.x; h];                         % append input from other neurons onto external stimulus

        % infer spike train for neuron 'i'
        fprintf('\nNeuron # %g\n',i)
        [S M II{i}.P]   = GOOPSI_main_v3_0(D(i).F,I{i}.P,Tim);
        II{i}.S.n       = S.n; 
        II{i}.S.h       = S.h; 
        II{i}.S.w_b     = S.w_b; 
        II{i}.M.nbar    = M.nbar;
    end

    % set inference for each neuron to the newly updated inference
    for i=1:Sim.Nc, I{i}.S = II{i}.S; I{i}.M = II{i}.M; end

    % given new inference for each neuron, update parameters
    for i=1:Sim.Nc

        % append external stimulus for neuron 'i' with spike histories from other cells
        h = zeros(Sim.Nc-1,Sim.T);                  % we append this to x to generate input into neuron from other neurons
        Pre=1:Sim.Nc;                               % generate list of presynaptic neurons
        Pre(Pre==i)=[];                             % remove self
        k=0;                                        % counter of dimension
        for j=Pre                                   % loop thru all presynaptic neurons
            k=k+1;                                  % generate input to neuron based on posterior mean spike train from neuron j
            h(k,:) = filter(1,[1 -(1-Sim.dt/P.tau_h)],I{j}.M.nbar);
        end
        Tim.x = [Sim.x; h];                         % append input from other neurons onto external stimulus

        fprintf('\nNeuron # %g\n',i)
        I{i}.P  = GOOPSI_Mstep_v1_0(Tim,I{i}.S,0,II{i}.P,D(i).F);
        EE{i}    = I{i}.P;
    end
    Phat{2}.omega = GetMatrix(Sim.Nc,EE);
    PlotMatrix(D,Sim,P,Phat)
    save(['SimTotal', num2str(Sim.T), '_', num2str(Sim.Nc), '_', num2str(tr)],'Phat','D','Sim','P')
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

% function PlotMatrix(D,Sim,P,Phat)
% 
% fig=figure(5); clf,
% col = [1 0 0; 0 0 1; 0 0.5 0; 1 0.5 0; 1 0 1];          % define colors for mean
% 
% tstart=2/Sim.dt;
% tdt=1/Sim.dt;
% tend = min(tstart+8*tdt,Sim.T);
% subplot(2,3,[1 2 3]), hold on
% for i=1:2
%     plot(z1(D(i).F(tstart:tend))+1,'Color',col(i,:));
% end
% axis('tight'),
% xticks  = tstart:tdt:tend;               % XTick positions
% set(gca,'YTick',[],'YTickLabel',[])
% set(gca,'XTick',xticks,'XTickLabel',xticks*Sim.dt)
% ylabel('Fluorescence')
% xlabel('Time (sec)')
% 
% clims(1)=min([P.omega(:)' Phat{1}.omega(:)' Phat{2}.omega(:)']);
% clims(2)=max([P.omega(:)' Phat{1}.omega(:)' Phat{2}.omega(:)']);
% subplot(234), imagesc(P.omega,clims), colormap(gray),
% set(gca,'XTick',[1:Sim.Nc],'YTick',[1:Sim.Nc]) %colorbar
% title('True matrix')
% ylabel('Presynaptic'), xlabel('Postsynaptic')
% 
% subplot(235), imagesc(Phat{1}.omega,clims), %colorbar
% set(gca,'XTick',[1:Sim.Nc],'YTick',[1:Sim.Nc]) %colorbar
% title('Matrix from spikes')
% % ylabel('Presynaptic'), xlabel('Postsynaptic')
% 
% subplot(236), imagesc(Phat{2}.omega,clims), %colorbar
% set(gca,'XTick',[1:Sim.Nc],'YTick',[1:Sim.Nc]) %colorbar
% title('Matrix from fluorescence')
% % ylabel('Presynaptic'), xlabel('Postsynaptic')
% 
% % print fig
% wh=[7 5];   %width and height
% set(fig,'PaperPosition',[0 11-wh(2) wh]);
% print('-depsc','SimConnector')
% print('-dpdf','SimConnector')
% 
% end