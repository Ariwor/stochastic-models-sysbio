% Generate sample paths using the tau-leap method with a fixed
% time-step tau for a network with four species and six reactions

function [X_results, jumps] = tau_leap(T_final)

% Initialize
tau = 0.01;
t=0;
X = [20 0 0 0];

X_results = [X];
jumps = [t];

c1=1;
c2=1;
c3=2;
c4=2;
c5=0.1;
c6=0.1;

% Define stoichiometry matrix
s1 = [-1 1 0 0];
s2 = [1 -1 0 0];
s3 = [-2 0 1 0];
s4 = [2 0 -1 0];
s5 = [1 0 0 -1];
s6 = [-1 0 0 1];
S_k = [s1; s2; s3; s4; s5; s6];

while t <= T_final

    % Calculate propensity functions
    a1 = c1*X(1);
    a2 = c2*X(2);
    a3 = c3*X(1)*((X(1) - 1)/2);
    a4 = c4*X(3);
    a5 = c5*X(1)*X(4);
    a6 = c6*X(1)*((X(1) - 1)/2);
    a_k = [a1 a2 a3 a4 a5 a6];

    % Draw samples
    N1 = poissrnd(a1*tau);
    N2 = poissrnd(a2*tau);
    N3 = poissrnd(a3*tau);
    N4 = poissrnd(a4*tau);
    N5 = poissrnd(a5*tau);
    N6 = poissrnd(a6*tau);
    N_k = [N1 N2 N3 N4 N5 N6];

    % Update time and number of molecules
    t = t + tau;
    for i=1:6
        X = X + S_k(i,:)*N_k(i);
    end

    % Check-Correct negative values
    check_1 = find(X<0);
    X(check_1(:))=0;

    X_results = [X_results; X];
    jumps=[jumps,t];
end
jumps=jumps(1:end-1);
X_results(end,:) = [];
end
