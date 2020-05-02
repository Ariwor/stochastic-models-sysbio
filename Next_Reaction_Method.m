% Generate sample paths using the Next Reaction Method
% for a network with four species and six reactions

% Based on the next reaction method described in:
% Anderson, D. F. (2007). A modified next reaction method for simulating chemical systems
% with time dependent propensities and delays. The Journal of chemical physics, 127(21), 214107.

function [X_results, jumps] = Next_Reaction_Method(T_final)

% Initialize
T_k = zeros(1,6);
P_k = zeros(1,6);

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

% Generate independent uniform random numbers
r1 = rand;
r2 = rand;
r3 = rand;
r4 = rand;
r5 = rand;
r6 = rand;

% Calculate P for every reaction
P1 = log(1/r1);
P2 = log(1/r2);
P3 = log(1/r3);
P4 = log(1/r4);
P5 = log(1/r5);
P6 = log(1/r6);

P_k = [P1 P2 P3 P4 P5 P6];

while t <= T_final
    
    % Calculate propensity functions
    a1 = c1*X(1);
    a2 = c2*X(2);
    a3 = c3*X(1)*((X(1) - 1)/2);
    a4 = c4*X(3);
    a5 = c5*X(1)*X(4);
    a6 = c6*X(1)*((X(1) - 1)/2);
    a_k = [a1 a2 a3 a4 a5 a6];

    % Calculate Delta_t
    Delta_t_1 = (P_k(1) - T_k(1))/a_k(1);
    Delta_t_2 = (P_k(2) - T_k(2))/a_k(2);
    Delta_t_3 = (P_k(3) - T_k(3))/a_k(3);
    Delta_t_4 = (P_k(4) - T_k(4))/a_k(4);
    Delta_t_5 = (P_k(5) - T_k(5))/a_k(5);
    Delta_t_6 = (P_k(6) - T_k(6))/a_k(6);
    Delta_t = [Delta_t_1, Delta_t_2, Delta_t_3, Delta_t_4 , Delta_t_5, Delta_t_6];

    % "Checking"
    check_inf = isinf(Delta_t);
    check_1 = find(check_inf==1);
    Delta_t(check_1(:))=10^12;

    Delta = min(Delta_t);
    mu = find(Delta_t==min(Delta_t));
    t = t + Delta;

    % Update time and number of molecules
    T_k = T_k + a_k*Delta;
    P_k(mu) = P_k(mu) + log(1/rand);
    X = X + S_k(mu,:);

    X_results = [X_results; X];
    jumps=[jumps,t];
end
jumps=jumps(1:end-1);
X_results(end,:) = [];
end