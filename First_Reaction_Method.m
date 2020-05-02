% Generate sample paths using the First Reaction Method
% for a network with four species and six reactions

function [X_results, jumps] = First_Reaction_Method(T_final)

% Initialize
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

    % Draw sample times
    T1 = -log(rand)/a1;
    T2 = -log(rand)/a2;
    T3 = -log(rand)/a3;
    T4 = -log(rand)/a4;
    T5 = -log(rand)/a5;
    T6 = -log(rand)/a6;
    T_k = [T1 T2 T3 T4 T5 T6];

    % "Checking"
    check_inf = isinf(T_k);
    check_1 = find(check_inf==1);
    T_k(check_1(:))=10^12;

    T_min = min(T_k);
    mu = find(T_k==min(T_k));

    % Update time and number of molecules
    t = t + T_min;
    X = X + S_k(mu,:);

    X_results = [X_results; X];
    jumps=[jumps,t];
end
jumps=jumps(1:end-1);
X_results(end,:) = [];
end
