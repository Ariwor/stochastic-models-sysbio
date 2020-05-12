% Generate sample paths using the Stochastic Simulation Algorithm (Gillespie’s Direct Method) for an example
% gene expression network, which incorporates transcriptional feedback from the protein molecules to enhance
% mRNA degradation. In other words, the transcription rate is equal to (k_r + g_fb*[protein number])*[mRNA number].

function [X_results, jumps] = SSA_Moment_Closure_ex(T_final)

% Initialize
t=0;
X = [0 0];
X_results = [X];
jumps = [t];

k_r=1;
k_p=5;
g_r=10;
g_p=2;
g_fb=6;

% Define stoichiometry matrix
s1 = [1 0];
s2 = [-1 0];
s3 = [0 1];
s4 = [0 -1];
S_k = [s1; s2; s3; s4];

while t <= T_final
    
% Calculate propensity functions
a1 = k_r;
a2 = (g_r + g_fb*X(end,2))*X(end,1);
a3 = k_p*X(end,1);
a4 = g_p*X(end,2);
a_k = [a1 a2 a3 a4];
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
% jumps=jumps(1:end-1);
% X_results(end,:) = [0 0];
end