function output = pop_oopsi(Cell,V0)


%% init nhat and Phat assuming each cell is independent

VV=V0;
Tmax            = V0.T; %min(500,length(Cell(1).F));
VV.T            = Tmax;
VV.Nspikehist   = 1;
VV.fast_iter_max= 0;
VV.smc_iter_max = 0;
VV.Ncells       = 1;
VV.Nparticles   = 50;

use_smc = 1;
cheat   = 1; % use true parameters
if use_smc
    VV.fast_do      = 0;
    VV.smc_do       = 1;
    for i=1:V0.Ncells
        VV.n    = Cell{i}.n;
        if cheat
            PP=Cell{i}.P;
            PP.omega= PP.omega(i,i);
            smc{i}  = run_oopsi(Cell{i}.F(1:Tmax),VV,PP);
        else
            smc{i}  = run_oopsi(Cell{i}.F(1:Tmax),VV);
        end
        Phat{i} = smc{i}.P;
        nhat{i} = smc{i}.E.nbar;
    end
else
    VV.smc_do   = 0;
    VV.fast_do  = 1;
    for i=1:V0.Ncells
        fast{i} = run_oopsi(Cell{i}.F(1:Tmax),VV);
        Phat{i} = fast{i}.P;
        nhat{i} = fast{i}.n;
    end
end

%% pop em pseudo-gibbs iteration

for i=1:V0.Ncells                               % for each cell

    fprintf('\nNeuron # %g\n',i)                % estimate weights

    % append external stimulus for neuron 'i' with spike histories from other cells
    h = zeros(V0.Ncells-1,V0.T);                % we append this to x to generate input into neuron from other neurons
    Pre=1:V0.Ncells;                            % generate list of presynaptic neurons
    Pre(Pre==i)=[];                             % remove self
    k=0;                                        % counter of dimension
    for j=Pre                                   % loop thru all presynaptic neurons
        k=k+1;                                  % generate input to neuron based on posterior mean spike train from neuron j
        if ~use_smc
            h(k,:) = filter(1,[1 -(1-V0.dt/Phat{i}.tau_h)],nhat{i});
        else
            h(k,:) = sum(smc{i}.E.h.*smc{i}.E.w);
        end
    end
    VV.x = [V0.x; [zeros(V0.Ncells-1,1) h(:,1:end-1)]];              % append input from other neurons onto external stimulus

    % set other variables for smc_oopsi_m_step
    VV.est_c = 0;
    VV.est_F = 0;
    VV.est_n = 1;
    VV.est_h = 1;
    VV.smc_plot=0;
    VV.StimDim = VV.Ncells;

    smc{i}.P.k = zeros(V0.Ncells+V0.StimDim-1,1);
    smc{i}.P.k(1:V0.StimDim) = Phat{i}.k;
    smc{i}.E.w_b=smc{i}.E.w;
    smc{i}.P  = smc_oopsi_m_step(VV,smc{i}.E,0,smc{i}.P,Cell{i}.F);
end

if use_smc
    output = smc;
else
    output = fast;
end

save('../../data/JoVo/sim1_inference','smc')

% if we have do the above recursively, then we will uncomment out this
% stuff, and put this whole section within an EM loop
%
% % for each neuron, infer spike train conditioned on previous EM
% % iterations spike history terms
% for i=1:V0.Ncells,                                 % infer spikes for each neuron
%
%     % append external stimulus for neuron 'i' with spike histories from other cells
%     h = zeros(V0.Ncells-1,V0.T);                  % we append this to x to generate input into neuron from other neurons
%     Pre=1:V0.Ncells;                               % generate list of presynaptic neurons
%     Pre(Pre==i)=[];                             % remove self
%     k=0;                                        % counter of dimension
%     for j=Pre                                   % loop thru all presynaptic neurons
%         k=k+1;                                  % generate input to neuron based on posterior mean spike train from neuron j
%         h(k,:) = filter(1,[1 -(1-V0.dt/Phat{i}.tau_h)],nhat{j});
%     end
%     VV.x = [V0.x; h];                         % append input from other neurons onto external stimulus
%
%     % infer spike train for neuron 'i'
%     fprintf('\nNeuron # %g\n',i)
%     if use_smc
%         [smc{i}] = run_oopsi(Cell{i}.F,VV,Phat{i});
% %         II{i}.S.n       = smc{i}.E.n;
% %         II{i}.S.h       = smc{i}.E.h;
% %         II{i}.S.w_b     = smc{i}.E.w;
% %         II{i}.M.nbar    = smc{i}.E.nbar;
%     end
% end

% set inference for each neuron to the newly updated inference
% for i=1:V0.Ncells, I{i}.S = II{i}.S; I{i}.M = II{i}.M; end

