clear all
fname='SimConnector40000_10_1.mat';
load(['/Users/joshyv/Research/projects/crcns_grant/' fname]);
D=R;

%%

Ts = [1000:1000:3000];
for j=1:length(Ts)

    % get estimate
    Sim.M       = 1;                                    % # spike history terms per neuron (fixed at one for this version of code)
    Sim.n_params= 1;                                    % if 1, estimate k
    Sim.h_params= 1;                                    % if 1, estimate omega (self-coupling)
    Sim.F_params= 0;                                    % if 1, estimate observation parameters
    Sim.C_params= 0;                                    % whether to compute
    Sim.StimDim = Sim.Nc;                               % set external stim dimesions to # cells
    Tim         = Sim;                                  % Tim is Sim for this estimation
    Tim.T       = Ts;                                   % number of spikes to use
    Tim.N       = 1;                                    % # of particles
    for i=1:Sim.Nc
        h=zeros(Sim.Nc-1,Sim.T);
        Pre=1:Sim.Nc;                                   % generate list of presynaptic neurons
        Pre(Pre==i)=[];                                 % remove self
        k=0;                                            % counter of dimension
        for j=Pre                                       % loop thru all presynaptic neurons
            k=k+1;                                      % generate input to neuron based on posterior mean spike train from neuron j
            h(k,:) = D(j).h;
        end
        Tim.x       = [Sim.x; h];                       % append input from other neurons onto external stimulus
        E           = P;
        E.omega     = E.omega(i,i);                     % initialize self-coupling term
        E.k         = E.k*ones(Sim.Nc,1);               % initialize external stim and cross-coupling terms
        D(i).w_b    = ones(1,Sim.T);
        Enew2{i}    = GOOPSI_Mstep_v1_0(Tim,D(i),0,E,D(i).F);
    end

    FiniteTime{j}.omega = GetMatrix(Sim.Nc,Enew2);
end




%%
nrows=5;
true_weights=P.omega(:);
n_weights=Phat{1}.omega(:);
OO=linspace(-2,1,2);
figure(4), clf,
for j=1:length(Ts)
    subplot(1,nrows,1), scatter(true_weights,n_weights)
    ylabel('true weights')
    xlabel('spike weights')
    hold on, plot(OO,OO,'k')
    n_xcorr=corrcoef(true_weights,n_weights);
    xx=0; for j=1:10, xx=xx+sum(R(j).n); end
    n_mean=xx/10;
    title(['<spikes/neuron>=' num2str(n_mean) ', r^2=' num2str(n_xcorr(2))])
end

subplot(1,nrows,j+1), scatter(true_weights,n_weights)
ylabel('true weights')
xlabel('spike weights')
hold on, plot(OO,OO,'k')
n_xcorr=corrcoef(true_weights,n_weights);
xx=0; for j=1:10, xx=xx+sum(R(j).n); end
n_mean=xx/10;
title(['<spikes/neuron>=' num2str(n_mean) ', r^2=' num2str(n_xcorr(2))])

    
F_weights=Phat{2}.omega(:);
subplot(1,nrows,nrows), scatter(true_weights,F_weights)
hold on, plot(OO,OO,'k')
ylabel('true weights')
xlabel('F weights')
F_xcorr=corrcoef(true_weights,F_weights);
title(['r^2=' num2str(F_xcorr(2))])

%     wh=[7 5];   %width and height
%     set(fnum,'PaperPosition',[0 11-wh(2) wh]);
print('-depsc',['scatter_' fname(13:end-4)])

