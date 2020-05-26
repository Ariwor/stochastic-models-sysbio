%%

clc;
clear all

% SSA

T_final = 10;
theta = 100;
P = 1000; % Number of paths

% Pre-allocation
X_results_T = zeros(P,1);

for k=1:P
[X_results, jumps] = SSA_autorepressive(T_final, theta);

% for T = 10
for index = 1:length(jumps)
    X_results_T(k,:) = X_results(end,:);
end
end

% Plotting
histogram(X_results_T, 'Normalization','pdf', 'Facecolor', '#7E2F8E')
title('Histogram of the number of X molecules at time = 10')
xlabel('X molecules')

% Saving
save('Data_T', 'X_results_T');

%%
clc;
clear all

theta_vec = [10: 200];

for i = theta_vec
    L(i) = FSP_MLE(i, 'Data_T.mat');
end

% Delete first 9 elements that were left empty

L(1:9)=[];

% Plotting

plot(theta_vec,L)
title('Plot of FSP-derived log-likehood function')
xlabel('theta')
ylabel('Log-likehood')
xlim([10 200])

% Find maximum L
% Note: We could also find the minimum if we used the negative FSP_MLE function

[L_max, index_max] = max(L);

% Find MLE for theta

theta_MLE = index_max + 9;

%% Alternative (general) solution
clc;
clear all

% Create handle for FSP_MLE function
FSP_MLE_handle = @(theta) -FSP_MLE(theta,'Data_T.mat');

% Optimization (without a specific range of values)
theta_MLE = fminsearch (FSP_MLE_handle, 100);

%%

% Parameters & Initialization

clc;
clear all

% Define length of the Markov chain trajectory
length = 10000;

% Set parameters
burn_in = 0.1*length; % Define burn-in (percentage of samples to discard)
freq = 4; % Define sub-sampling frequency
samples_stored = (0.9*length)/(freq + 1); % Define remaining trajectory to sub-sample and store

% Initial (starting) condition (point)
theta = 100; 

% Storing
Theta = zeros(samples_stored,1); % samples to store from the remaining trajectory 

% Mean & variance of the proposal (lognormal) distribution
mu = 0;
sigma = 1;

%% Metropolis-Hastings algorithm

% Burn-in run

for i = 1:burn_in
    [theta] = MH_FSP(theta, sigma); 
end

% Sub-sampling and storing with the defined frequency

for i = 1:samples_stored
    for j = 1:freq
        [theta] = MH_FSP(theta, sigma); 
    end
    Theta(i) = theta; % store every 5-th sample
end

%% Plotting

histogram(Theta, 'Normalization','pdf', 'Facecolor', '#7E2F8E')
title('Histogram (posterior distribution)')
xlabel('theta')
ylabel('Propability')