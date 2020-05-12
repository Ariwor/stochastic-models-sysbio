% Approximate, using stochastic differential equations, a concentration process X_omega = X(t)/omega = (mRNA(t)/omega, Protein(t)/omega),
% where omega is the volume of the system, for an example gene expression network which incorporates transcriptional feedback from the protein
% molecules to enhance mRNA degradation.

% Two approximations are employed, the Diffusion or Langevin Approximation and the Van Kampen or Linear Noise Approximation.
% Both of them make use of Brownian motion.

% Number of trajectories
P = 5; 

omega = 1000;
k_r = 1;
g_r = 10;
k_p = 5;
g_p = 2;
g_fb = 6;

% Initialize
t=0;
Z = [0 0];
Z_results = [Z];
time = [t];

% Define stoichiometry matrix
s1 = [1 0];
s2 = [-1 0];
s3 = [0 1];
s4 = [0 -1];
S_k = [s1; s2; s3; s4]';

% Define timestep
Delta_t = 0.1;

%%

% Solve the stochastic differential equation (SDE) for the Diffusion
% (Langevin) Approximation of the concentration process X_omega, using the
% Euler-Maruyama scheme

for j=1:P
    
for i=0:99
    
    t_i = i*Delta_t;
    
    % Calculate propensity functions for w_k(z(t))
    w_k_1 = k_r;
    w_k_2 = (g_r + g_fb*Z(end, 2))*Z(end,1);
    w_k_3 = k_p*Z(end,1);
    w_k_4 = g_p*Z(end,2);
    w_k = [w_k_1; w_k_2; w_k_3; w_k_4];

    % Draw independent  normal  random  variables  with  mean  0  and variance 1
    h_1 = normrnd(0,1);
    h_2 = normrnd(0,1); 
    h_3 = normrnd(0,1);
    h_4 = normrnd(0,1);
    h_k = [h_1; h_2; h_3; h_4]';
    
    % Euler-Maruyama scheme
    Z = Z + ((w_k(1)*S_k(:,1) + w_k(2)*S_k(:,2) + w_k(3)*S_k(:,3) + w_k(4)*S_k(:,4))*Delta_t)' + (sqrt((Delta_t)/omega)*(S_k(:,1)*sqrt(w_k(1))*h_k(1) + S_k(:,2)*sqrt(w_k(2))*h_k(2) + S_k(:,3)*sqrt(w_k(3))*h_k(3) + S_k(:,4)*sqrt(w_k(4))*h_k(4)))';
    
    Z_results = [Z_results; Z];
    time = [time, t_i];
    i=i+1;
end
cum_Z_results_mRNA(:,j) = Z_results(:,1);
cum_Z_results_Protein(:,j) = Z_results(:,2);
Z_results = [0 0];

% Plotting diffusion approximation for the mRNA species
plot(time, cum_Z_results_mRNA(:,j));
title('Trajectories (N=5) of the diffusion approximation Z for the mRNA species')
xlabel('t') 
hold on
legend({'Trajectory 1','Trajectory 2','Trajectory 3', 'Trajectory 4', 'Trajectory 5'},'Location','southwest')
time_stored = time;
time = [0];
end

%%
hold off
% Plotting the diffusion approximation for the Protein species
for j=1:P
plot(time_stored, cum_Z_results_Protein(:,j));
hold on
title('Trajectories (N=5) of the diffusion approximation Z for the Protein species')
xlabel('t') 
legend({'Trajectory 1','Trajectory 2','Trajectory 3', 'Trajectory 4', 'Trajectory 5'},'Location','southwest')
end

%%

% Solve the ordinary differenential equations for the deterministic
% concentration process (x(t))(for large volume omega, X_omega converges to this
% deterministic process).

% Integrating the ODEs
tspan = [0 10];
x0 = [0 0];

[t,x] = ode15s(@ODEs_deterministic_process, tspan, x0);

% Plotting 

tiledlayout(2,1)

nexttile
plot(t,x(:,1));
title('Trajectory of the mRNA species')
xlabel('t') 
legend({'mRNA copy-number','Location','southwest'})
nexttile
plot(t,x(:,2));
title('Trajectory of the Protein species')
xlabel('t') 
legend({'Protein copy-number','Location','southwest'})

%%

% Solve the differential equations for the deterministic (x(t)) and fluctuation (L(t)) process using the Euler-Maruyama scheme
% to estimate the Linear Noise Approximation (LNA) for the concentration process X_omega

% Initialize

t=0;
x = [0 0];
L = [0 0];
LNA = [0 0];

x_results = [x];
L_results = [L];
LNA_results = [LNA];
time = [t];

for j=1:P
for i=0:99
    t_i = i*Delta_t;
    
    % Calculate propensity functions for w_k(x(t))
    w_k_1 = k_r;
    w_k_2 = (g_r + g_fb*x(end,2))*x(end,1);
    w_k_3 = k_p*x(end,1);
    w_k_4 = g_p*x(end,2);
    w_k = [w_k_1; w_k_2; w_k_3; w_k_4];

    % Draw independent  normal  random  variables  with  mean  0  and variance 1
    h_1 = normrnd(0,1);
    h_2 = normrnd(0,1); 
    h_3 = normrnd(0,1);
    h_4 = normrnd(0,1);
    h_k = [h_1; h_2; h_3; h_4]';
    
    % Euler-Maruyama scheme
    x = x + ((w_k(1)*S_k(:,1) + w_k(2)*S_k(:,2) + w_k(3)*S_k(:,3) + w_k(4)*S_k(:,4))*Delta_t)';
    L = L + ([-(g_r + g_fb*x(end,2)), -g_fb*x(end,1); k_p, - g_p]*(L)'*Delta_t)' + (sqrt(Delta_t)*(S_k(:,1)*w_k(1)*h_k(1) + S_k(:,2)*w_k(2)*h_k(2) + S_k(:,3)*w_k(3)*h_k(3) + S_k(:,4)*w_k(4)*h_k(4)))';
    
    % Calculate LNA
    LNA = x + L/sqrt(omega);
    
    x_results = [x_results; x];
    L_results = [L_results; L];
    LNA_results = [LNA_results; LNA];
    time = [time, t_i];
    i=i+1;
end

cum_LNA_results_mRNA(:,j) = LNA_results(:,1);
cum_LNA_results_Protein(:,j) = LNA_results(:,2);

LNA_results = [0 0];

% Plotting LNA approximation of the mRNA species
plot(time, cum_LNA_results_mRNA(:,j));
title('Trajectories (N=5) of the LNA approximation for the mRNA species')
xlabel('t') 
hold on
legend({'Trajectory 1','Trajectory 2','Trajectory 3', 'Trajectory 4', 'Trajectory 5'},'Location','southwest')
time_stored = time;
time = [0];

end

%%
hold off
% Plotting LNA approximation of the Protein species
for j=1:P
plot(time_stored, cum_LNA_results_Protein(:,j));
hold on
title('Trajectories (N=5) of the LNA approximation for the Protein species')
xlabel('t') 
legend({'Trajectory 1','Trajectory 2','Trajectory 3', 'Trajectory 4', 'Trajectory 5'},'Location','southwest')
end

%%

% Find the covariance matrix C for the fluctuation process L(t) by solving
% the respective Lyapunov equation

x_eq = [0.0883, 0.2208]; % Steady-state values of X(t)
A = [-(g_r + g_fb*x_eq(2)), -g_fb*x_eq(1); k_p, -g_p];
B = [sqrt(k_r), -sqrt((g_r + g_fb*x_eq(2))*x_eq(1)), 0, 0; 0, 0, sqrt(k_p*x_eq(1)), -sqrt(g_p*x_eq(2))];
Q = B*B';

% Solve Lyapunov equation and find the covariance matrix C
C = lyap(A,Q);
