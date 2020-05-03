% Estimate mean, variance and coefficient of determination for each species
% at a specific time, after simulation with the First Reaction Method, SSA (Direct Method), Next
% Reaction Method and tau-leaping with fixed tau

%% Number of paths
P = 10000; 

%% First Reaction Method
for k=1:P
[X_results, jumps] = First_Reaction_Method(5);
% for T=5
for w_FRM = 1:length(jumps)
    K(k,:) = X_results(end,:);
end
end
mean_FRM = mean(K);
var_FRM = var(K);
CV_FRM = std(K)./mean(K);

%% Stochastic Simulation Algorithm (Gillespie’s Direct Method)
for k=1:P
[X_results, jumps] = SSA(5);
% for T=5
for w_SSA = 1:length(jumps)
    L(k,:) = X_results(end,:);
end
end
mean_SSA = mean(L);
var_SSA = var(L);
CV_SSA = std(L)./mean(L);

%% Next Reaction Method
for k=1:P
[X_results, jumps] = Next_Reaction_Method(5);
% for T=5
for w_NRM = 1:length(jumps)
    M(k,:) = X_results(end,:);
end
end
mean_NRM = mean(M);
var_NRM = var(M);
CV_NRM = std(M)./mean(M);

%% tau-leap with fixed tau
for k=1:P
[X_results, jumps] = tau_leap(5);
% for T=5
for w_tau = 1:length(jumps)
    Z(k,:) = X_results(end,:);
end
end
mean_tau = mean(Z);
var_tau = var(Z);
CV_tau = std(Z)./mean(Z);
