function S = Pop_forward_v1_0(Sim,F,P)
% the function does the forward part of the pop pf
% here, we are simply doing vanilla particle filtering, ie, sample using
% the prior.
%
% Inputs---
% Sim:  simulation metadata
% F:    fluorescence (Nc x T)
% P:    current parameter estimates
%
% Outputs---
% S:    simulation states
%%

% extize particle info
Z  = zeros(Sim.Np,Sim.T);                  % matrix of zeros that is Np x T
S(1).p  = Z;                   % extize rate
S(1).n  = Z;                   % extize spike counts
S(1).C  = P.C_init*(1+Z);           % extize calcium
S(1).h  = Z;                  % extize spike history terms
S(1).eps_c = sqrt(P.sig2_c)*randn(Sim.Np,Sim.T);% generate noise on C
S(1).eps_h = sqrt(P.sig2_h)*randn(Sim.Np,Sim.T);% generate noise on C
S(1).U_sampl = rand(Sim.Np,Sim.T);                % random samples
for i=1:Sim.Nc, S(i)=S(1); end

% preprocess stuff for setting weights and stratified resampling
w_f   = 1/Sim.Np*(1+Z);            % extize forward weights
Neff  = Sim.Np*(1+Z(1,:));              % extize N_{eff}
ints        = linspace(0,1,Sim.Np+1);            % generate intervals
diffs       = ints(2)-ints(1);                  % generate interval size
U_resamp    = repmat(ints(1:end-1),Sim.T,1)+diffs*rand(Sim.T,Sim.Np); % resampling matrix

% extize misc stuff
col = [1 0 0; 0 0 1];          % define colors for mean

% do the particle filter
for t=2:Sim.T
    for i=1:Sim.Nc

        % generate samples
        S(i).h(:,t) = P.g*S(i).h(:,t-1)+S(i).n(:,t-1) + S(i).eps_h(:,t);
        y_t         = P.kx(t);
        for j=1:Sim.Nc
            y_t     = y_t + P.omega(i,j)*S(j).h(:,t);
        end
        S(i).p(:,t) = 1-exp(-exp(y_t)*Sim.dt);            % update rate for those particles with y_t<0
        S(i).n(:,t) = S(i).U_sampl(:,t)<S(i).p(:,t);      % sample n
        S(i).C(:,t) = (1-P.a)*S(i).C(:,t-1)+P.A*S(i).n(:,t)+P.a*P.C_0+S(i).eps_c(:,t);% sample C

        % compute weights for each cell
        F_mu(i,:)   = P.alpha*Hill_v1(P,S(i).C(:,t))+P.beta;   % compute E[F_t]
        F_var(i,:)  = P.gamma*Hill_v1(P,S(i).C(:,t))+P.zeta;   % compute V[F_t]
        ln_w(i,:)   = -0.5*(F(i,t)-F_mu(i,:)).^2./F_var(i,:) - 0.5*log(F_var(i,:));% compute log of weights
    end

    % compute weights for each particle
    ln_wsum = sum(ln_w);
    ln_wsum = ln_wsum+log(w_f(:,t-1)');   % update log(weights)
    ln_wsum = ln_wsum-max(ln_wsum);                       % subtract the max to avoid rounding errors
    w       = exp(ln_wsum);                            % exponentiate to get actual weights
    w_f(:,t)= w/sum(w);                             % normalize to define a legitimate distribution

    % stratified resample
    Neff(t)  = 1/sum(w_f(:,t).^2);    % store N_{eff}
    if Neff(t) < Sim.Np/2                         % if weights are degenerate or we are doing prior sampling then resample
        [foo,ind]   = histc(U_resamp(t,:),[0  cumsum(w_f(:,t))']);
        [ri,ri]     = sort(rand(Sim.Np,1));      % these 3 lines stratified resample
        ind         = ind(ri);

        for i=1:Sim.Nc
            S(i).p(:,t)   = S(i).p(ind,t);      % resample probabilities (necessary?)
            S(i).n(:,t)   = S(i).n(ind,t);      % resample calcium
            S(i).C(:,t)   = S(i).C(ind,t);      % resample calcium
            S(i).h(:,t)   = S(i).h(ind,t);      % resample h's
        end
        w_f(:,t) = 1/Sim.Np*(1+Z(:,1)); % reset weights
    end

    % plot stuff
    figure(2), clf
    for i=1:Sim.Nc
        subplot(511), plot(z1(F(i,1:t))+1,'Color',col(i,:)), hold on, axis('tight'), ylabel('F')
        subplot(512), stem(mean(S(i).n(:,1:t)),'Color',col(i,:)), hold on, axis('tight'), ylabel('n')
        subplot(513), plot(mean(S(i).p(:,2:t)),'Color',col(i,:)), hold on, axis('tight'), ylabel('p')
    end
    subplot(514), plot(w_f(:,1:t)'), axis('tight'), ylabel('w')
    subplot(515), plot(Neff(1:t)), ylabel('Neff')

end %for time loop