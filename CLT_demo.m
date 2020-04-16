% This code demonstrates the Central Limit Theorem. X is a continuous random variable
% with logistic distribution. The aim is to show that the distribution of Z
% (defined in the code) approaches the standard normal distribution, with increasing n.

%% Compute one realization of Z with a desired n value

n = 100;
X = zeros(1, 100);
for i=1:n
    u = rand(1); % generate u from the uniform distribution
    X(i) = 2 - log((1-u)/u)*sqrt(3); % inverse function for an example c.d.f. (of the logistic distribution)
    Z = (sum(X) - 2*n)/(pi*sqrt(n));
end
%% Compute m realizations of Z with n = 100

n = 100;
m = 1000;
X = zeros(1,100);
Z = zeros(1,1000);
for j=1:m
for i=1:n
    u = rand(1); % generate u from the uniform distribution
    X(i) = 2 - log((1-u)/u)*sqrt(3); % inverse function for an example c.d.f. (of the logistic distribution)
end
    Z(j) = (sum(X) - 2*n)/(pi*sqrt(n));
end

% Plotting
edges = [-10,10];
h1 = histogram(Z, edges, 'Normalization','probability', 'FaceColor','#77AC30');
h1.BinWidth = 1;
title('Histogram of random variable Z');

%% Compute m realizations of Z with n = [100, 1000, 10000]

n_v = [100, 1000, 10000];
m = 1000;
Z = zeros(1,1000);
for n = n_v
for j=1:m
    for i=1:n
        u = rand(1); % generate u from the uniform distribution
        X(i) = 2 - log((1-u)/u)*sqrt(3); % inverse function for an example c.d.f. (of the logistic distribution)
    end
    Z(j) = (sum(X) - 2*n)/(pi*sqrt(n));
end

% Plotting
figure();
tiledlayout(2,1)
nexttile
edges = [-10,10];
h1 = histogram(Z, edges, 'Normalization','probability', 'FaceColor','#77AC30');
h1.BinWidth = 1;

% Plotting the p.d.f. of a standard normal distribution
x = [-10:.1:10];
pdf_n = normpdf(x);
nexttile
plot(x, pdf_n, 'r');
end
