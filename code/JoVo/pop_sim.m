% clear, clc

%% 1) set simulation metadata

V.T       = 2000;               % # of time steps
V.dt      = 1/60;               % time step size
V.StimDim = 1;                  % # dimensions of external stimulus
V.x       = ones(V.StimDim,V.T);% stimulus
V.Ncells  = 3;                  % # of cells
V.plot    = 0;

%% 2) initialize parameters

rate        = 5;                            % expected spike rate (assuming no spike history terms and V.x=1)
P.k         = log(-log(1-rate*V.dt)/V.dt);  % linear filter
P.k         = P.k*ones(V.StimDim,1);        % initialize k to the right number of dimensions
P.tau_c     = 0.5;                          % calcium decay time constant (sec)
P.A         = 20;                           % jumps size (\mu M)
P.C_0       = 20;                           % baseline [Ca++] (\mu M)
P.C_init    = P.C_0;                          % initial [Ca++] (\mu M)
P.sigma_c   = 1;                            % noise on
P.n         = 1.0;                          % hill equation exponent
P.k_d       = 200;                          % hill coefficient
P.alpha     = 1;                            % F_max
P.beta      = 0;                            % F_min
P.gamma     = 1e-6;                         % scaled variance
P.zeta      = 4*P.gamma;                    % constant variance
P.tau_h     = 0.05;                         % time constant
P.sigma_h   = 0.01;                         % stan dev of noise

% make ring network weights
w=1;
P.omega=diag(-2*ones(V.Ncells,1));
P.omega(1,2)=w;
P.omega(end,end-1)=w/2;
for i=2:2:V.Ncells-1
    P.omega(i,i-1)=w/2;
    P.omega(i,i+1)=w;
end
for i=3:2:V.Ncells-1
    P.omega(i,i-1)=-w;
    P.omega(i,i+1)=-w*2;
end

for i=1:V.Ncells; Cell{i}.P=P; end

%% 3) simulate data

Cell{1}.h     = zeros(1,V.T);
Cell{1}.n     = zeros(1,V.T);
Cell{1}.C     = zeros(1,V.T);
Cell{1}.F     = zeros(1,V.T);
for ii=2:V.Ncells; Cell{ii}=Cell{1}; end

kx      = P.k'*V.x;                                     % external input to neuron
eps_c   = P.sigma_c*sqrt(V.dt)*randn(V.Ncells,V.T);     % generate noise on calcium
U_sampl = rand(V.Ncells,V.T);                           % generate random number to use for sampling
eps_h   = repmat(P.sigma_h*sqrt(V.dt),V.Ncells,V.T).*randn(V.Ncells,V.T); % generate noise on spike history
eps_F   = randn(V.Ncells,V.T);                          % generate noise on fluorescence
p       = zeros(V.Ncells,V.T);                          % prob of spiking for each cell at each tim
y       = zeros(V.Ncells,V.T);                          % input to each cell at each time

for t=2:V.T                                             % update states
    for i=1:V.Ncells                                    % loop over presynaptic cells
        Cell{i}.h(t)   = (1-V.dt./P.tau_h).*Cell{i}.h(t-1) + Cell{i}.n(t-1) + eps_h(i,t); % update h terms
    end

    for i=1:V.Ncells                                    % loop over presynaptic cells
        y(i,t)      = kx(t);                            % initialize input to cell
        for j=1:V.Ncells                                % loop of post-synaptic cells
            y(i,t)  = y(i,t)+P.omega(i,j)*Cell{j}.h(t); % generate operand for rate function
        end
        p(i,t)          = 1-exp(-exp(y(i,t))*V.dt);     % generate prob of spiking
        Cell{i}.n(t)    = U_sampl(i,t)<p(i,t);         % sample from bernoulli with prob p_t
        Cell{i}.C(t)    = (1-V.dt/P.tau_c)*Cell{i}.C(t-1) + (V.dt/P.tau_c)*P.C_0 + P.A*Cell{i}.n(t) + eps_c(i,t);   % update calcium
        s               = Hill_v1(P,Cell{i}.C(t));     % compute saturated calcium
        Cell{i}.F(t)    = (P.alpha*s+P.beta)+sqrt(P.gamma*s+P.zeta).*eps_F(i,t); % update fluorescence
        if Cell{i}.F(t) <= 0; Cell{i}.F(t) = eps; end % keep fluorescence positive
    end
end

%% 4) plot simulation results

if V.plot==1
    col = jet(V.Ncells);          % define colors for mean
    figure(1), clf
    for i=1:V.Ncells
        h1=subplot(311); hold on, plot(((Cell{i}.F)./max(Cell{i}.F))+1,'Color',col(i,:)); stem(Cell{i}.n,'Color',col(i,:)), axis('tight'), ylabel('F')
        h2=subplot(312); hold on, plot(y(i,:),'Color',col(i,:)), axis('tight'), ylabel('y') %stem(Cell{i}.n,'Color',col(i,:)),
        h3=subplot(313); hold on, plot(p(i,:),'Color',col(i,:)), axis('tight'), ylabel('p')
    end
    linkaxes([h1 h2 h3],'x')
    legend('1','2')
%     for i=1:V.Ncells, disp(sum(Cell{i}.n)); end
end

save('../../data/sim1_data','Cell','V')