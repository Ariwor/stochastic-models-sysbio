% This code demonstrates a biological example/phenomenon that gives rise to the binomial distribution.
% A mother cell (with m plasmids) divides, giving rise to two daughter cells. X is a random variable
% that is defined as the number of plasmids that are inherited by the first daughter cell (with probability p)
% after division. It is demonstrated that X follows a binomial distribution (m,p).

%% Set up number of plasmids m and probability p at time of division

m = 10;
p = 0.6;

% Generate random variables u, uniformly distributed over [0,1]
u = rand(1,m);

% Calculate number of plasmids X
X = sum(u<p);

%% Compute n = 10000 realizations of X
tic
n = 10000;
X = zeros(1,10000);
for i=1:n
    u = rand(1,m);
    X(i) = sum(u<p);
end

% Plotting
tiledlayout(2,1)
nexttile
histogram(X,'Normalization','probability', 'FaceColor','#77AC30');
title('Histogram of random variable X');

% Plotting the p.m.f of the Binomial distribution
x = [0:10];
pmf_b = binopdf(x,m,p);
nexttile
scatter(x,pmf_b,'r');
title('Probability mass function of the Binomial distribution');
toc