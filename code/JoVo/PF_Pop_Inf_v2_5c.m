function PF_Pop_Inf_v2_5c(T,pl)

%% simulate or load data

pop_sim

use_spikes=1;
if use_spikes
    Phat{1}.omega = GetWSpikes(Cell,V,P);
end
%% set up structures

V.Nspikehist=1;
V.smc_plot=0;
V.smc_iter_max=1;
V.Nparticles=50;
V.fast_do=0;
V.smc_do=1;

V.est_n=1;
V.est_h=1;
V.est_F=0;
V.est_c=0;

P.lik       = 0;
Tim         = V;                        % Tim is V for this estimation
Tim.N       = 1;                        % # of particles
Tim         = V;                        % copy V structure for input to Mstep function
E           = P;
E.omega     = E.omega(1,1);             % initialize self-coupling term
E.k         = E.k*ones(V.StimDim,1);    % initialize external stim and cross-coupling terms

for i=1:V.Ncells,
    I{i}.M.nbar = zeros(1,V.T);         % initialize spike trains
    I{i}.P      = E;                    % initialize parameters
    smc{i}.P    = E;
end

%% infer spikes and parameters with independent assumption

for i=1:V.Ncells                    % infer spikes for each neuron

    % append external stimulus for neuron 'i' with spike histories from other cells
    h = zeros(V.Ncells-1,V.T);      % we append this to x to generate input into neuron from other neurons
    Pre=1:V.Ncells;                 % generate list of presynaptic neurons
    Pre(Pre==i)=[];                 % remove self
    k=0;                            % counter of dimension
    for j=Pre                       % loop thru all presynaptic neurons
        k=k+1;                      % generate input to neuron based on posterior mean spike train from neuron j
        h(k,:) = filter(1,[1 -(1-V.dt/P.tau_h)],I{j}.M.nbar);
    end
    Tim.x = [V.x; h];               % append input from other neurons onto external stimulus

    % infer spike train for neuron 'i'
    fprintf('\nNeuron # %g\n',i)
    smc{i}          = run_oopsi(Cell{i}.F,V,smc{i}.P);
    II{i}.S.n       = smc{i}.E.n;
    II{i}.S.h       = smc{i}.E.h;
    II{i}.S.w_b     = smc{i}.E.w;
    II{i}.M.nbar    = smc{i}.E.nbar;
    II{i}.P         = smc{i}.P;
    smc{i}.E.w_b    = smc{i}.E.w;
end

% set inference for each neuron to the newly updated inference
for i=1:V.Ncells,
    I{i}.S = II{i}.S;
    I{i}.M = II{i}.M;
end

%% infer connectivity given inferred spikes

for i=1:V.Ncells

    % append external stimulus for neuron 'i' with spike histories from other cells
    h = zeros(V.Ncells-1,V.T);      % we append this to x to generate input into neuron from other neurons
    Pre=1:V.Ncells;                 % generate list of presynaptic neurons
    Pre(Pre==i)=[];                 % remove self
    k=0;                            % counter of dimension
    for j=Pre                       % loop thru all presynaptic neurons
        k=k+1;                      % generate input to neuron based on posterior mean spike train from neuron j
        h(k,:) = filter(1,[1 -(1-V.dt/P.tau_h)],I{j}.M.nbar);
    end
    Tim.x = [V.x; h];               % append input from other neurons onto external stimulus

    fprintf('\nNeuron # %g\n',i)
    II{i}.P.k=smc{i}.P.k*ones(V.Ncells,1);
    smc{i}.P = smc_oopsi_m_step(Tim,smc{i}.E,0,II{i}.P,Cell{i}.F);
    EE{i}    = smc{i}.P;
end
Phat{2}.omega = GetMatrix(V.Ncells,EE);
PlotMatrix(Cell,V,P,Phat)
save(['../../data/VConnector', num2str(V.T), '_', num2str(V.Ncells)],'Phat','Cell','V','P','smc','II','I','E')
Fs=1024; ts=0:1/Fs:1; sound(sin(2*pi*ts*200)),

end



function PlotMatrix(Cell,V,P,Phat)

fig=figure(5); clf,
col = [1 0 0; 0 0 1; 0 0.5 0; 1 0.5 0; 1 0 1];          % define colors for mean

tstart=2/V.dt;
tdt=1/V.dt;
tend = min(tstart+8*tdt,V.T);
subplot(2,3,[1 2 3]), hold on
for i=1:2
    plot(((Cell{i}.F(tstart:tend))./max(Cell{i}.F(tstart:tend)))+i,'Color',col(i,:));
end
axis('tight'),
xticks  = tstart:tdt:tend;               % XTick positions
set(gca,'YTick',[],'YTickLabel',[])
set(gca,'XTick',xticks,'XTickLabel',xticks*V.dt)
ylabel('Fluorescence')
xlabel('Time (sec)')

clims(1)=min([P.omega(:)' Phat{1}.omega(:)' Phat{2}.omega(:)']);
clims(2)=max([P.omega(:)' Phat{1}.omega(:)' Phat{2}.omega(:)']);
% clims(1)=min([P.omega(:)' Phat{2}.omega(:)']); %Phat{2}.omega(:)']);
% clims(2)=max([P.omega(:)' Phat{2}.omega(:)']);% Phat{2}.omega(:)']);
subplot(234), imagesc(P.omega,clims), colormap(gray),
set(gca,'XTick',[1:V.Ncells],'YTick',[1:V.Ncells]) %colorbar
title('True matrix')
ylabel('Presynaptic'), xlabel('Postsynaptic')

subplot(235),
imagesc(Phat{1}.omega,clims), %colorbar
set(gca,'XTick',[1:V.Ncells],'YTick',[1:V.Ncells]) %colorbar
title('Matrix from spikes')
% ylabel('Presynaptic'), xlabel('Postsynaptic')

subplot(236), imagesc(Phat{2}.omega,clims), %colorbar
set(gca,'XTick',[1:V.Ncells],'YTick',[1:V.Ncells]) %colorbar
title('Matrix from fluorescence')
% ylabel('Presynaptic'), xlabel('Postsynaptic')

% print fig
wh=[7 5];   %width and height
set(fig,'PaperPosition',[0 11-wh(2) wh]);
print('-depsc','VConnector')
print('-dpdf','VConnector')

end