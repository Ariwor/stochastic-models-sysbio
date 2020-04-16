% Simulation of a simple mRNA transcription/degradation model (birth-death process)

%% First Reaction Method simulation - Option 1

T = 10;
lambda = 2;
gamma = 1;
P = 5; % Number of paths
t = zeros(1,100);
X = zeros(1,100);
figure();
hold on
for k=1:P
    j = 1;  
    t(1)= -log(rand)./lambda;
    X(1)=1;
    while t(j) <= T
        gamma_X = gamma*X(j);
        time_birth = -log(rand)./lambda;
        if X(j) == 0
            time_death = 10^6;
        else
            time_death = -log(rand)./gamma_X;
        end
            if time_birth < time_death
                X(j+1) = X(j) + 1;
                t(j+1) = t(j) + time_birth;
            else
                X(j+1) = X(j) - 1;
                t(j+1) = t(j) + time_death;
            end
            j=j+1;
    end
    stairs(t(1:j-1), X(1:j-1));
end

%% First Reaction Method simulation - Option 2

T = 10;
lambda = 2;
gamma = 1;
P=5; % Number of paths
mRNA = zeros(1,100);
jumps = zeros(1,100);
figure();
hold on
for k=1:P
    X = 0;
    t = 0;
    jumps = [0];
    mRNA = [X];
    while t <= T
        t = t - log(rand) / (lambda + X*gamma);
        if t > T
            break;
        else
            jumps = [jumps t];
            if rand < lambda / (lambda + X*gamma)
                X = X + 1;
            else
                X = X - 1;
            end
        mRNA = [mRNA X];
        end
    end
    stairs(jumps, mRNA);
end

%% First Reaction Method simulation and plot realizations of mRNA populations at specific times

T = 10;
lambda = 2;
gamma = 1;
P = 10000; % Number of paths
t = zeros(1,100);
X = zeros(1,100);
K = zeros(1,10000);
L = zeros(1,10000);
M = zeros(1,10000);
for k=1:P
    j = 1;  
    t(1)= -log(rand)./lambda;
    X(1)=1;
    while t(j) <= T
        gamma_X = gamma*X(j);
        time_birth = -log(rand)./lambda;
        if X(j) == 0
            time_death = 10^6;
        else
            time_death = -log(rand)./gamma_X;
        end
            if time_birth < time_death
                X(j+1) = X(j) + 1;
                t(j+1) = t(j) + time_birth;
            else
                X(j+1) = X(j) - 1;
                t(j+1) = t(j) + time_death;
            end
            j=j+1;
    end
% for T=2
    for w_2 = 1:length(t) 
        y_2 = find(t>2);
        if y_2(1) == 1
            K(k) = X(y_2(1));
        else
            K(k) = X(y_2(1)-1);
        end
    end
% for T=5
    for w_5 = 1:length(t)
        y_5 = find(t>5);
        if y_5(1) == 1
            L(k) = X(y_5(1));
        else
            L(k) = X(y_5(1)-1);
        end
    end
% for T=10
    for w_10 = 1:length(t)
        y_10 = find(t>10);
        if y_10(1) == 1
            M(k) = X(y_10(1));
        else
            M(k) = X(y_10(1)-1);
        end
    end
end

% Plotting
figure();
tiledlayout(3,1)
nexttile
h1 = histogram(K, 'Normalization','probability', 'FaceColor','#77AC30');
nexttile
h2 = histogram(L, 'Normalization','probability', 'FaceColor','#77AC30');
nexttile
h3 = histogram(M, 'Normalization','probability', 'FaceColor','#77AC30');

%% Estimation of mean, variance and coefficient of variation

% For T=2
mean_2 = mean(K);
var_2 = var(K);
CV_2 = std(K)/mean(K);

% For T=5
mean_5 = mean(L);
var_5 = var(L);
CV_5 = std(L)/mean(L);

% For T=10
mean_10 = mean(M);
var_10 = var(M);
CV_10 = std(M)/mean(M);