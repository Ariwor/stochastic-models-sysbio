% Generate sample paths of the birth-death process using the random time-change (RTC)
% representation of a birth-death network

% Based on the next reaction method described in:
% Anderson, D. F. (2007). A modified next reaction method for simulating chemical systems
% with time dependent propensities and delays. The Journal of chemical physics, 127(21), 214107.

function [X, jumps] = RTC_next_reaction_method_birth_death(kappa, gamma, T_final, X_0)

% Initialize
T1 = 0;
T2 = 0;

t=0;
X = [X_0];
jumps = [t];

% Calculate propensity functions, for each reaction
a1 = kappa;
a2 = gamma*X(end);

% Generate independent uniform random numbers
r1 = rand;
r2 = rand;

% Calculate P, the first firing time of Y, for every reaction
P1 = log(1/r1);
P2 = log(1/r2);

while t <= T_final

    % Calculate Delta_t
    Delta_t_1 = (P1 - T1)/a1;
    Delta_t_2 = (P2 - T2)/a2;

    % Update time and number of molecules
    if Delta_t_1 < Delta_t_2
        % Birth reaction fires
        X = [X, X(end) + 1];
        t = t + Delta_t_1;
        P1 = P1 + log(1/rand);
        T1 = T1 + a1*Delta_t_1;
        T2 = T2 + a2*Delta_t_1;

    else
        % Death reaction fires
        X = [X, X(end) - 1];
        t = t + Delta_t_2;
        P2 = P2 + log(1/rand);
        T1 = T1 + a1*Delta_t_2;
        T2 = T2 + a2*Delta_t_2;

    end
    a2 = gamma*X(end);
    jumps=[jumps,t];
    
end
X=X(1:end-1);
jumps=jumps(1:end-1);
end