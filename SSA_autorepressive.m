% Generate sample paths using the Stochastic Simulation Algorithm (Gillespie’s Direct Method
% for an example auto-repressive gene-expression network.

function [X_results, jumps] = SSA_autorepressive(T_final, theta)

% Initialize
t=0;
X = [0];
X_results = [X];
jumps = [t];

alpha = 2;
gamma = 1;

% Define stoichiometry matrix
s1 = +1;
s2 = -1;
S_k = [s1; s2];

while t <= T_final
    
% Calculate propensity functions

lamda_x = theta/(1+(X(end,1))^(alpha));
birth = lamda_x;
death = gamma*X(end,1);
a_k = [birth death];
a0 = sum(a_k);

% Draw sample time for next reaction

t0 = -log(rand)/a0;

% Draw sample identity for next reaction

a_cumsum = cumsum(a_k)/a0;
i0 = min(find(rand <= a_cumsum, 1, 'first'));

% Update time and number of molecules
t = t + t0;
X = X + S_k(i0,:);

X_results = [X_results; X];
jumps=[jumps,t];

end
end