% Generate sample paths using the Stochastic Simulation Algorithm (Gillespie’s Direct Method) for the Schlogl model.
% Rate constants depend on a volume parameter V, which is defined by the user.

function [X_results, jumps] = SSA_Schlog_vol(T_final, V)

% Initialize
t=0;
X = [0 0];
X_results = [X];
jumps = [t];

c_1 = V;
c_2 = 1;
c_3 = 5/V;
c_4 = 0.2/V;
c_5 = 5;

% Define stoichiometry matrix
s1 = [1 0];
s2 = [2 -1];
s3 = [-1 +1];
s4 = [-1 0];
s5 = [-1 0];
S_k = [s1; s2; s3; s4; s5];

while t <= T_final
    
% Calculate propensity functions
a1 = c_1;
a2 = c_2*X(end,2);
a3 = c_3*(X(end,1)*(X(end,1)-1))/2;
a4 = c_4*X(end,1)*X(end,2);
a5 = c_5*X(end,1);
a_k = [a1 a2 a3 a4 a5];
a0 = sum(a_k);

% Draw sample time for next reaction
t0 = -log(rand)/a0;

% Draw sample identity for next reaction
a_cumsum = cumsum(a_k)/a0;
i0 = min(find(rand <= a_cumsum));

% Update time and number of molecules
t = t + t0;
X = X + S_k(i0,:);

X_results = [X_results; X];
jumps=[jumps,t];
end
%jumps=jumps(1:end-1);
%X_results(end,:) = [0 0];
end