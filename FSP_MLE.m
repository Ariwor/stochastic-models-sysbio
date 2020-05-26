% A function to compute log-likelihood for the data observations 'Data' given a
% parameter 'theta', using the finite state projection (FSP) algorithm for an example
% auto-repressive gene expression network. The user inputs are a data file and a value of 'theta'.

function L = FSP_MLE(theta, Data)

% Define the time point of the data samples
T = 10;

% Import the data samples at a specific time point
D = importdata(Data);

% Define the propensity functions
alpha = 2;
gamma = 1;
birth = @(X,theta) theta/(1+(X)^(alpha));
death = @(X,gamma) gamma*X;

% Define the truncated state-space
trunc_space_upper = 100;
trunc_space_max = trunc_space_upper+1;

% Pre-allocation
A = sparse(zeros(trunc_space_max,trunc_space_max));

% Create the FSP matrices

% Matrix A
for i = 1:trunc_space_max
    
    % Calculate propensities
    index = i-1;
    w1 = birth(index, theta);
    w2 = death(index, gamma);
    
    % Fill the diagonal entries of A
    A (i,i) = -w1-w2;
    
    % Fill the off diagonal entries of A
    if (index<trunc_space_upper)
        A (i+1, i) = w1;
    end
    
    if (index>0)
        A (i-1, i) = w2;
    end
end

% Matrix B
B = -sum(A);

% Matrix C
C = [A zeros(numel(A(:,1)),1);B 0];

% Define initial distribution (conditions)
q0 = [1;zeros(trunc_space_max,1)];

% Define initial  distribution (uniform between 0 and 9)
% q0 = [0.1*ones(10,1);zeros(trunc_space_max-9,1)];

% Solve the system at time T
P = expm(C*T)*q0;

% Compute the log-likehood
L = sum(log(P(D+1)));

end