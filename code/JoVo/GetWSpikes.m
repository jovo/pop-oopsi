function omega = GetWSpikes(Cell,V,P)
% estimate connection matrix directly from spikes

V.est_n=1;
V.est_h=1;
V.est_F=0;
V.est_c=0;
V.Nspikehist=1;
V.Nparticles=1;
V.smc_plot=0;
V.smc_iter_max=1;


P.lik       = 0;
Tim         = V;                                  % Tim is V for this estimation
Tim.N       = 1;                                    % # of particles
for i=1:V.Ncells
    h=zeros(V.Ncells-1,V.T);
    Pre=1:V.Ncells;                                   % generate list of presynaptic neurons
    Pre(Pre==i)=[];                                 % remove self
    k=0;                                            % counter of dimension
    for j=Pre                                       % loop thru all presynaptic neurons
        k=k+1;                                      % generate input to neuron based on posterior mean spike train from neuron j
        h(k,:) = Cell{j}.h;
    end
    Tim.x       = [V.x; h];                       % append input from other neurons onto external stimulus
    E           = P;
    E.omega     = E.omega(i,i);                     % initialize self-coupling term
    E.k         = E.k*ones(V.Ncells,1);               % initialize external stim and cross-coupling terms
    Cell{i}.w_b = ones(1,V.T);
    Enew2{i}    = smc_oopsi_m_step(Tim,Cell{i},0,E,Cell{i}.F);
end

omega = GetMatrix(V.Ncells,Enew2);

% if pl==1
figure(3), clf,
clims(1)=min(min(P.omega(:)),min(omega(:)));
clims(2)=max(max(P.omega(:)),max(omega(:)));
subplot(121), imagesc(P.omega,clims), colormap(gray), %colorbar
subplot(122), imagesc(omega,clims), %colorbar
%     Fs=1024; ts=0:1/Fs:1; sound(sin(2*pi*ts*200)),
%keyboard
% end